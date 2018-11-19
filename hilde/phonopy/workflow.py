""" Provide a full highlevel phonopy workflow """

from os.path import isfile
from pathlib import Path
import pickle
import numpy as np

from ase.calculators.socketio import SocketIOCalculator

from hilde.helpers.converters import dict2atoms
from hilde.helpers.paths import cwd
from hilde.phonon_db.row import PhononRow
import hilde.phonopy.phono as ph
from hilde.trajectory.phonopy import metadata2file, step2file
from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.trajectory import reader as traj_reader
from hilde.helpers.k_grid import update_k_grid

def initialize_phonopy_attach_calc(atoms, calc, supercell_matrix, displacement=0.01, workdir="."):
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)
    for sc in scs:
        sc.calc = calc
    return phonon, supercell, scs

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
    force_constants_file="force_constants.dat",
    supercell_file="geometry.in.supercell",
    pickle_file="phonon.pick",
    #    fingerprint_file="fingerprint.dat",
):
    """ perform a full phonopy calculation """

    trajectory = Path(trajectory).absolute()

    watchdog = Watchdog(walltime=walltime)

    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)

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

    force_sets = []
    with cwd(workdir, mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:
            # save inputs and supercell
            metadata2file(atoms, calc, phonon, trajectory)
            supercell.write(supercell_file, format="aims", scaled=True)

            if socketio_port is not None:
                calc = iocalc

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

    postprocess(
        phonon,
        trajectory=trajectory,
        workdir=workdir,
        force_constants_file=force_constants_file,
        pickle_file=pickle_file,
        #        supercell_file=supercell_file,
        #        fingerprint_file=fingerprint_file,
    )


def postprocess(
    phonon,
    calculated_atoms=None,
    trajectory="phonopy_trajectory.yaml",
    workdir=".",
    force_constants_file="force_constants.dat",
    #    supercell_file="geometry.in.supercell",
    #    fingerprint_file="fingerprint.dat",
    displacement=0.01,
    fireworks=False,
    pickle_file="phonon.pick"
    #    **kwargs,
):
    """ Phonopy postprocess """

    if fireworks:
        phonon = PhononRow(phonon).to_phonon()
        phonon.generate_displacements(
            distance=displacement, is_plusminus="auto", is_diagonal=True
        )

    if calculated_atoms:
        if fireworks:
            calculated_atoms = [dict2atoms(cell) for cell in calculated_atoms]
        calculated_atoms = sorted(
            calculated_atoms, key=lambda x: x.info["displacement_id"]
        )
    elif isfile(trajectory):
        calculated_atoms = traj_reader(trajectory)
    else:
        raise ValueError("Either calculated_atoms or trajectory must be defined")

    force_sets = [atoms.get_forces() for atoms in calculated_atoms]

    # compute and save force constants
    force_constants = ph.get_force_constants(phonon, force_sets)
    np.savetxt(Path(workdir) / force_constants_file, force_constants)

    with open(pickle_file, "wb") as fp:
        pickle.dump(phonon, fp)
