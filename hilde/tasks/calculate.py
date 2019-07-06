"""
 Functions to run several related calculations with using trajectories as cache
"""
from pathlib import Path
import numpy as np
from ase.calculators.socketio import SocketIOCalculator
from hilde.helpers import talk
from hilde.helpers.compression import backup_folder as backup
from hilde.helpers.socketio import get_port
from hilde.helpers.watchdogs import SlurmWatchdog as Watchdog
from hilde.trajectory import (
    get_hashes_from_trajectory,
    hash_atoms,
    metadata2file,
    step2file,
)

from hilde.son import son

from hilde.helpers.paths import cwd
from hilde.helpers.aims import get_aims_uuid_dict

calc_dirname = "calculations"


def cells_and_workdirs(cells, base_dir):
    """generate tuples of atoms object and workingdirectory path

    Parameters
    ----------
    cells: list of ase.atoms.Atoms
        The cells to assign workdirs to
    base_dir: str
        Path to the base working directory

    Yields
    ------
    cell: ase.atoms.Atoms
        The particular cell
    workdir: Path
        The working directory for cell
    """
    for ii, cell in enumerate(cells):
        workdir = Path(base_dir) / f"{ii:05d}"
        yield cell, workdir


def calculate(atoms, calculator, workdir="."):
    """Perform a dft calculation with ASE

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The structure to calculate
    calculator: ase.calculators.calulator.Calculator:
        The calculator to used to get the properties
    workdir: str or Path
        Path to the working directory

    Returns
    -------
    calc_atoms: ase.atoms.Atoms
        atoms with all properties calculated
    """

    if atoms is None:
        return atoms
    calc_atoms = atoms.copy()
    calc_atoms.calc = calculator
    with cwd(workdir, mkdir=True):
        calc_atoms.calc.calculate(calc_atoms)
    return calc_atoms


def calculate_socket(
    atoms_to_calculate,
    calculator,
    metadata=None,
    settings=None,
    trajectory="trajectory.son",
    workdir="calculations",
    save_input=False,
    backup_folder="backups",
    backup_after_calculation=True,
    check_settings_before_resume=True,
    **kwargs,
):
    """perform calculations for a set of atoms objects, while able to use the socket

    Parameters
    ----------
    atoms_to_calculate: list of ase.atoms.Atoms
        list with atoms to calculate
    calculator: ase.calculators.calulator.Calculator
        calculator to use
    metadata: dict
        metadata information to store to trajectory
    settings: dict
        the settings used to set up the calculation
    trajectory: str or Path
        path to write trajectory to
    workdir: str or Path
        working directory
    backup_folder: str or Path
        directory to back up calculations to
    check_settings_before_resume: bool
        only resume when settings didn't change

    Returns
    -------
    bool
        Wether all structures were computed or not

    Raises
    ------
    RuntimeError
        If the lattice changes significantly
    """

    # create watchdog
    watchdog = Watchdog(buffer=1)

    # create working directories
    workdir = Path(workdir).absolute()
    trajectory = workdir / trajectory
    backup_folder = workdir / backup_folder
    calc_dir = workdir / calc_dirname

    # perform backup if calculation folder exists
    backup(calc_dir, target_folder=f"{backup_folder}")

    # save fist atoms object for computation
    atoms = atoms_to_calculate[0].copy()

    # handle the socketio
    socketio_port = get_port(calculator)
    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calculator

    # append settings to metadata
    if settings:
        metadata["settings"] = settings.to_dict()
        if save_input:
            with cwd(workdir, mkdir=True):
                geometry_file = Path(settings.geometry.file)
                print(geometry_file)
                if not geometry_file.exists():
                    settings.atoms.write(str(geometry_file), format="aims")
                settings.obj["workdir"] = workdir
                settings.write()

    # fetch list of hashes from trajectory
    precomputed_hashes = get_hashes_from_trajectory(trajectory)

    # perform calculation
    n_cell = -1
    with cwd(calc_dir, mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:

            # log metadata and sanity check
            if check_settings_before_resume:
                try:
                    old_metadata, _ = son.load(trajectory)
                    check_metadata(metadata, old_metadata)
                    talk(f"resume from {trajectory}")
                except FileNotFoundError:
                    metadata2file(metadata, trajectory)

            if socketio_port is not None:
                calc = iocalc
            else:
                calc = calculator

                # fix for EMT calculator
                try:
                    calc.initialize(atoms)
                    talk("calculator initialized.")
                except AttributeError:
                    pass

            for n_cell, cell in enumerate(atoms_to_calculate):
                # skip if cell is None or already computed
                if cell is None:
                    continue
                if hash_atoms(cell) in precomputed_hashes:
                    talk(f"Structure {n_cell + 1} already computed, skip.")
                    continue

                # perform backup of settings
                if settings:
                    settings.write()

                # make sure a new calculation is started
                calc.results = {}
                atoms.calc = calc

                # update calculation_atoms and compute force
                atoms.info = cell.info
                atoms.positions = cell.positions

                # check the cell
                try:
                    diff = np.linalg.norm(atoms.cell.array - cell.cell.array)
                except AttributeError:
                    diff = np.linalg.norm(atoms.cell - cell.cell)
                if diff > 1e-10:
                    msg = "lattice has changed by {diff}, FIXME!"
                    raise RuntimeError(msg)

                talk(f"Compute structure {n_cell + 1} of {len(atoms_to_calculate)}")
                atoms.calc.calculate(atoms, system_changes=["positions"])

                # peek into aims file and grep for uuid
                meta = get_aims_uuid_dict()

                step2file(atoms, atoms.calc, trajectory, metadata=meta)

                if watchdog():
                    break

    # backup
    if backup_after_calculation:
        backup(calc_dir, target_folder=f"{backup_folder}")

    if n_cell < len(atoms_to_calculate) - 1:
        return False

    return True


def check_metadata(new_metadata, old_metadata, keys=["calculator"]):
    """check if metadata sets coincide and sanity check geometry

    Parameters
    ----------
    new_metadata: dict
        The metadata of the current calculation
    old_metadata: dict
        The metadata from the trajectory.son file
    keys: list of str
        Keys to check if the metadata agree with

    Raises
    ------
    ValueError
        If the keys do not coincide
    """
    nm = new_metadata
    om = old_metadata

    if "atoms" in new_metadata:
        new_atoms = new_metadata["atoms"]
        old_atoms = old_metadata["atoms"]
        if new_atoms["symbols"] != old_atoms["symbols"]:
            msg = f"Compare Structures:\n{new_atoms}\n{old_atoms}"
            raise RuntimeError(msg)

    for key in keys:
        if key is "walltime":
            continue
        if isinstance(nm[key], dict):
            check_metadata(nm[key], om[key], keys=nm[key].keys())
        if nm[key] != om[key]:
            msg = f"{key} changed: from {nm[key]} to {om[key]}"
            raise ValueError(msg)
