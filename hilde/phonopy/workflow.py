""" Provide a full highlevel phonopy workflow """

from pathlib import Path

from ase.calculators.socketio import SocketIOCalculator

from hilde.settings import Settings
from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.paths import cwd, move_to_dir
import hilde.phonopy.wrapper as ph
from hilde.trajectory.phonons import metadata2dict, step2file, last_from_yaml, to_yaml
from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.compression import backup_folder as backup
from hilde.helpers.socketio import get_port
from .postprocess import postprocess


_calc_dirname = "calculations"


def last_calculation_id(trajectory):
    """ return the id of the last computed supercell """
    disp_id = -1

    try:
        dct = last_from_yaml(trajectory)
        disp_id = dct["phonopy"]["id"]
    except (FileNotFoundError, KeyError):
        pass

    return disp_id


def run(
    atoms,
    calc,
    supercell_matrix,
    kpt_density=None,
    displacement=0.01,
    symprec=1e-5,
    trajectory="trajectory.yaml",
    walltime=1800,
    workdir=".",
    primitive_file="geometry.in.primitive",
    supercell_file="geometry.in.supercell",
    force_constants_file="force_constants.dat",
    pickle_file="phonon.pick",
    backup_folder="backups",
    db_path=None,
    **kwargs,
):

    calculation_atoms, calc, supercell, scs, phonon, metadata = preprocess(
        atoms, calc, supercell_matrix, kpt_density, displacement
    )

    calculate(
        atoms=calculation_atoms,
        calculator=calc,
        primitive=atoms,
        supercell=supercell,
        supercells_with_displacements=scs,
        metadata=metadata,
        trajectory=trajectory,
        walltime=walltime,
        workdir=workdir,
        primitive_file=primitive_file,
        supercell_file=supercell_file,
    )

    postprocess(
        phonon,
        trajectory=trajectory,
        workdir=workdir,
        force_constants_file=force_constants_file,
        pickle_file=pickle_file,
        db_path=db_path,
        **kwargs,
    )

    return True


def preprocess(
    atoms, calc, supercell_matrix, kpt_density=None, displacement=0.01, symprec=1e-5
):

    # Phonopy preprocess
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement, symprec)

    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    # grab first geometry for computation
    calculation_atoms = scs[0].copy()

    if kpt_density is not None:
        update_k_grid(calculation_atoms, calc, kpt_density)

    metadata = metadata2dict(atoms, calc, phonon)

    return calculation_atoms, calc, supercell, scs, phonon, metadata


def calculate(
    atoms,
    calculator,
    primitive,
    supercell,
    supercells_with_displacements,
    metadata,
    trajectory="trajectory.yaml",
    walltime=1800,
    workdir=".",
    primitive_file="geometry.in.primitive",
    supercell_file="geometry.in.supercell",
    backup_folder="backups",
    **kwargs,
):
    """ perform a full phonopy calculation """
    print(workdir)
    workdir = Path(workdir)
    trajectory = (workdir / trajectory).absolute()
    backup_folder = workdir / backup_folder
    calc_dir = workdir / _calc_dirname

    # take the literal settings for running the task
    settings = Settings()

    watchdog = Watchdog(walltime=walltime, buffer=1)

    # handle the socketio
    socketio_port = get_port(calculator)
    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calculator

    # save input geometries and settings
    with cwd(workdir, mkdir=True):
        primitive.write(primitive_file, format="aims", scaled=True)
        supercell.write(supercell_file, format="aims", scaled=False)
        settings.write()

    # perform calculation
    n_cell = -1
    with cwd(calc_dir, mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:
            if last_calculation_id(trajectory) < 0:
                to_yaml(metadata, trajectory, mode="w")

            if socketio_port is not None:
                calc = iocalc
            else:
                calc = calculator

            for n_cell, cell in enumerate(supercells_with_displacements):
                if "forces" in calc.results:
                    del(calc.results['forces'])
                atoms.calc = calc

                if cell is None:
                    continue

                # if precomputed:
                if n_cell <= last_calculation_id(trajectory):
                    continue

                # update calculation_atoms and compute force
                atoms.positions = cell.positions
                _ = atoms.get_forces()

                step2file(atoms, atoms.calc, n_cell, trajectory)

                if watchdog():
                    break

    # backup and restart if necessary
    if n_cell < len(supercells_with_displacements) - 1:
        backup(calc_dir, target_folder=backup_folder)
        return False

    return True


def preprocess_fireworks(atoms, calc, supercell_matrix, displacement=0.01):
    """ phonopy preprocess returning supercells with attached calculator for FW """
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True
    for sc in scs:
        sc.calc = calc
    return phonon, supercell, scs, calc
