from pathlib import Path
import numpy as np

from ase.calculators.socketio import SocketIOCalculator

from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.paths import cwd
from hilde.trajectory.phonopy import metadata2file, step2file
import hilde.phonopy.phono as ph


def phonopy(
    atoms,
    calc,
    supercell_matrix,
    displacement=0.01,
    trajectory="phonopy_trajectory.yaml",
    socketio_port=None,
    walltime=1800,
    workdir=".",
    force_constants_file="force_constants.dat",
    supercell_file="geometry.in.supercell",
    fingerprint_file="fingerprint.dat",
):
    """ perform a full phonopy calculation """

    watchdog = Watchdog(walltime=walltime)

    trajectory = Path(trajectory).absolute()

    if "compute_forces" in calc.parameters:
        calc.parameters["compute_forces"] = True

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc

    # get a correctly shaped supercell matrix
    if np.size(supercell_matrix) == 1:
        supercell_matrix = supercell_matrix * np.eye(3)
    elif np.size(supercell_matrix) == 3:
        supercell_matrix = np.diag(supercell_matrix)
    elif np.size(supercell_matrix) == 9:
        supercell_matrix = np.asarray(supercell_matrix).reshape((3, 3))
    else:
        raise Exception(
            f"Supercell matrix must have 1, 3, 9 elements, has {len(supercell_matrix)}"
        )

    # preprocess
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)

    # grab first geometry for computation
    calculation_atoms = scs[0].copy()

    force_sets = []

    with cwd(workdir, mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:
            if socketio_port is not None:
                calc = iocalc

            calculation_atoms.calc = calc

            # save inputs and supercell
            metadata2file(atoms, calc, phonon, trajectory)
            supercell.write(supercell_file, format="aims", scaled=True)

            for ii, cell in enumerate(scs):

                if cell is None:
                    continue

                # update calculation_atoms and compute force
                calculation_atoms.positions = cell.positions
                force_sets.append(calculation_atoms.get_forces())

                step2file(calculation_atoms, calculation_atoms.calc, ii, trajectory)

                if watchdog():
                    break

        # compute and save force constants
        force_constants = ph.get_force_constants(phonon, force_sets)
        np.savetxt(force_constants_file, force_constants)
