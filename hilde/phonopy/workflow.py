""" Provide a full highlevel phonopy workflow """

from pathlib import Path

from ase.calculators.socketio import SocketIOCalculator

from hilde.settings import Settings
from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.paths import cwd, move_to_dir
import hilde.phonopy.wrapper as ph
from hilde.trajectory.phonopy import metadata2file, step2file, last_from_yaml
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


def phonopy(
    atoms,
    calc,
    supercell_matrix,
    kpt_density=None,
    displacement=0.01,
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
    #    fingerprint_file="fingerprint.dat",
):
    """ perform a full phonopy calculation """

    workdir = Path(workdir)
    trajectory = (workdir / trajectory).absolute()
    backup_folder = workdir / backup_folder
    calc_dir = workdir / _calc_dirname

    # take the literal settings for running the task
    settings = Settings()

    # make sure forces are computed
    if calc.name == 'aims':
        calc.parameters['compute_forces'] = True

    watchdog = Watchdog(walltime=walltime, buffer=1)

    # Phonopy preprocess
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)

    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    socketio_port = get_port(calc)
    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc

    # grab first geometry for computation
    calculation_atoms = scs[0].copy()

    if kpt_density is not None:
        update_k_grid(calculation_atoms, calc, kpt_density)

    # save input geometries and settings
    with cwd(workdir, mkdir=True):
        atoms.write(primitive_file, format="aims", scaled=True)
        supercell.write(supercell_file, format="aims", scaled=False)
        settings.write()

    # perform calculation
    n_cell = -1
    with cwd(calc_dir, mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:
            if last_calculation_id(trajectory) < 0:
                metadata2file(atoms, calc, phonon, trajectory)

            if socketio_port is not None:
                calc = iocalc

            for n_cell, cell in enumerate(scs):
                calculation_atoms.calc = calc

                if cell is None:
                    continue

                # if precomputed:
                if n_cell <= last_calculation_id(trajectory):
                    continue

                # update calculation_atoms and compute force
                calculation_atoms.positions = cell.positions
                _ = calculation_atoms.get_forces()

                step2file(calculation_atoms, calculation_atoms.calc, n_cell, trajectory)

                if watchdog():
                    break

    # backup and restart if necessary
    if n_cell < len(scs) - 1:
        backup(calc_dir, target_folder=backup_folder)
        return False

    postprocess(
        phonon,
        trajectory=trajectory,
        workdir=workdir,
        force_constants_file=force_constants_file,
        pickle_file=pickle_file,
        db_path=db_path,
        **kwargs,
        # fingerprint_file=fingerprint_file,
    )

    return True


def initialize_phonopy_attach_calc(atoms, calc, supercell_matrix, displacement=0.01):
    """ phonopy preprocess returning supercells with attached calculator for FW """
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)
    if calc.name == 'aims':
        calc.parameters['compute_forces'] = True
    for sc in scs:
        sc.calc = calc
    return phonon, supercell, scs, calc
