from os.path import isfile
from pathlib import Path
import numpy as np

from ase.calculators.socketio import SocketIOCalculator

from hilde.helpers.converters import dict2results, dict2atoms
from hilde.helpers.paths import cwd
from hilde.phonon_db.row import PhononRow
import hilde.phonopy.phono as ph
from hilde.trajectory.phonopy import metadata2file, step2file
from hilde.watchdogs import WallTimeWatchdog as Watchdog

def initialize_phonopy(
    atoms,
    calc,
    supercell_matrix=None,
    displacement=0.01,
    socketio_port=None,
    **kwargs
):
    """ perform a full phonopy calculation """
    if supercell_matrix is None:
        raise ValueError("The supercell matrix must be defined")

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
    return ph.preprocess(atoms, supercell_matrix, displacement)


def analyze_phonpy(
    phonon,
    calculated_atoms=None,
    trajectory="phonopy_trajectory.yaml",
    workdir=".",
    force_constants_file="force_constants.dat",
    supercell_file="geometry.in.supercell",
    fingerprint_file="fingerprint.dat",
    displacement=0.01,
    fireworks=False,
    **kwargs
):
    if fireworks:
        phonon = PhononRow(phonon).to_phonon()
        phonon.generate_displacements(
            distance=displacement, is_plusminus="auto", is_diagonal=True
        )
    if calculated_atoms:
        if fireworks:
            calculated_atoms = [dict2atoms(cell) for cell in calculated_atoms]
        calculated_atoms = sorted(calculated_atoms, key=lambda x: x.info['displacement_id'])
    elif isfile(trajectory):
        traj = from_yaml(trajectory)
        calculated_atoms = [dict2results(el) for el in traj[1:]]
    else:
        raise ValueError("Either calculated_atoms or trajectory must be defined")
    force_sets = [atoms.get_forces() for atoms in calculated_atoms]
    with cwd(workdir):
        # compute and save force constants
        force_constants = ph.get_force_constants(phonon, force_sets)
        np.savetxt(force_constants_file, force_constants)

def phonopy(
    atoms,
    calc,
    supercell_matrix=None,
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
    phonon, supercell, scs = initialize_phonopy(atoms, calc, supercell_matrix, displacement, socketio_port)

    # grab first geometry for computation
    calculation_atoms = scs[0].copy()

    force_sets = []

    with cwd(workdir, mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:
            if socketio_port is not None:
                calc = iocalc

            # save inputs and supercell
            metadata2file(atoms, calc, phonon, trajectory)
            supercell.write(supercell_file, format="aims", scaled=True)

            for ii, cell in enumerate(scs):
                calculation_atoms.calc = calc

                if cell is None:
                    continue

                # update calculation_atoms and compute force
                calculation_atoms.positions = cell.positions
                force_sets.append(calculation_atoms.get_forces())

                step2file(calculation_atoms, calculation_atoms.calc, ii, trajectory)

                if watchdog():
                    break
        # # compute and save force constants
        # force_constants = ph.get_force_constants(phonon, force_sets)
        # np.savetxt(force_constants_file, force_constants)

    analyze_phonpy(
        phonon,
        trajectory,
        workdir,
        force_constants_file,
        supercell_file,
        fingerprint_file
    )
