""" Provide a full highlevel phonopy workflow """

import sys
from pathlib import Path

from ase.calculators.socketio import SocketIOCalculator

from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.paths import cwd
import hilde.phonopy.wrapper as ph
from hilde.trajectory.phonopy import metadata2file, step2file, last_from_yaml
from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.compression import backup_folder
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
    trajectory="phonopy_trajectory.yaml",
    socketio_port=None,
    walltime=1800,
    workdir=".",
    primitive_file="geometry.in.primitive",
    supercell_file="geometry.in.supercell",
    force_constants_file="force_constants.dat",
    bandstructure_file="bandstructure.pdf",
    pickle_file="phonon.pick",
    db_path=None,
    backup_settings=None,
    restart_cmd=None,
    **kwargs,
    #    fingerprint_file="fingerprint.dat",
):
    """ perform a full phonopy calculation """

    trajectory = Path(trajectory).absolute()
    calc_dir = Path(workdir) / _calc_dirname

    watchdog = Watchdog(walltime=walltime, buffer=1)

    if Path(workdir).exists():
        print(f"***Caution, directory {workdir} exists, chance of data loss.")

    # Phonopy preprocess
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)

    # make sure forces are computed (aims only)
    if "compute_forces" in calc.parameters:
        calc.parameters["compute_forces"] = True

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc

    # grab first geometry for computation
    calculation_atoms = scs[0].copy()

    if kpt_density is not None:
        update_k_grid(calculation_atoms, calc, kpt_density)

    # save input geometries and configuration
    with cwd(workdir, mkdir=True):
        atoms.write(primitive_file, format="aims", scaled=True)
        supercell.write(supercell_file, format="aims", scaled=True)

        if backup_settings is not None:
            backup_settings.write()

    # perform calculation
    with cwd(calc_dir, mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:
            # save inputs and supercell
            if last_calculation_id(trajectory) < 0:
                metadata2file(atoms, calc, phonon, trajectory)

            if socketio_port is not None:
                calc = iocalc

            for ii, cell in enumerate(scs):
                calculation_atoms.calc = calc

                if cell is None:
                    continue

                # if precomputed:
                if ii <= last_calculation_id(trajectory):
                    continue

                # update calculation_atoms and compute force
                calculation_atoms.positions = cell.positions
                _ = calculation_atoms.get_forces()

                step2file(calculation_atoms, calculation_atoms.calc, ii, trajectory)

                if watchdog():
                    break

    # backup and restart if necessary
    if ii < len(scs) - 1:
        backup_folder(calc_dir, target_folder=workdir)
        if restart_cmd is not None:
            print(f"restart with {restart_cmd}")
            run(restart_cmd.split())
        else:
            sys.exit("**Watchdog: running out of time! No restart script.")

    postprocess(
        phonon,
        trajectory=trajectory,
        workdir=workdir,
        force_constants_file=force_constants_file,
        bandstructure_file=bandstructure_file,
        pickle_file=pickle_file,
        db_path=db_path,
        **kwargs,
        # fingerprint_file=fingerprint_file,
    )


def initialize_phonopy_attach_calc(atoms, calc, supercell_matrix, displacement=0.01):
    """ phonopy preprocess returning supercells with attached calculator for FW """
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)
    for sc in scs:
        sc.calc = calc
    return phonon, supercell, scs
