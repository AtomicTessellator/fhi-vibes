""" Provide a full highlevel phonopy workflow """

from pathlib import Path
import pickle
import numpy as np

from ase.calculators.socketio import SocketIOCalculator

from hilde.helpers.converters import dict2atoms
from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.paths import cwd
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import PhononRow
import hilde.phonopy.wrapper as ph
from hilde.phonopy import displacement_id_str
from hilde.structure.convert import to_Atoms
from hilde.trajectory.phonopy import metadata2file, step2file
from hilde.trajectory import reader as traj_reader
from hilde.watchdogs import WallTimeWatchdog as Watchdog


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
    **kwargs,
    #    fingerprint_file="fingerprint.dat",
):
    """ perform a full phonopy calculation """

    trajectory = Path(trajectory).absolute()
    workdir = Path(workdir).absolute()

    watchdog = Watchdog(walltime=walltime, buffer=1)

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

    # save input geometries
    with cwd(workdir, mkdir=True):
        atoms.write(primitive_file, format="aims", scaled=True)
        supercell.write(supercell_file, format="aims", scaled=True)

    # perform calculation
    with cwd(workdir / "calculations", mkdir=True):
        with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc:
            # save inputs and supercell
            metadata2file(atoms, calc, phonon, trajectory)

            if socketio_port is not None:
                calc = iocalc

            for ii, cell in enumerate(scs):
                calculation_atoms.calc = calc

                if cell is None:
                    continue

                # update calculation_atoms and compute force
                calculation_atoms.positions = cell.positions
                _ = calculation_atoms.get_forces()

                step2file(calculation_atoms, calculation_atoms.calc, ii, trajectory)

                if watchdog():
                    raise ResourceWarning("**Watchdog: running out of time!")

    postprocess(
        phonon,
        trajectory=trajectory,
        workdir=workdir,
        force_constants_file=force_constants_file,
        bandstructure_file=bandstructure_file,
        pickle_file=pickle_file,
        **kwargs,
        # fingerprint_file=fingerprint_file,
    )


def postprocess(
    phonon,
    calculated_atoms=None,
    trajectory="phonopy_trajectory.yaml",
    workdir=".",
    force_constants_file="force_constants.dat",
    bandstructure_file="bandstructure.pdf",
    displacement=0.01,
    fireworks=False,
    pickle_file="phonon.pick",
    db_path=None,
    **kwargs,
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
            calculated_atoms, key=lambda x: x.info[displacement_id_str]
        )
    elif Path(trajectory).is_file():
        calculated_atoms = traj_reader(trajectory)
    else:
        raise ValueError("Either calculated_atoms or trajectory must be defined")

    force_sets = [atoms.get_forces() for atoms in calculated_atoms]

    # compute and save force constants
    force_constants = ph.get_force_constants(phonon, force_sets)
    np.savetxt(Path(workdir) / force_constants_file, force_constants)

    with (Path(workdir) / pickle_file).open("wb") as fp:
        pickle.dump(phonon, fp)
    if db_path:
        update_phonon_db(
            db_path,
            to_Atoms(phonon.get_unitcell()),
            phonon,
            calc_type="phonons",
            symprec=phonon._symprec,
            sc_matrix_2=list(phonon.get_supercell_matrix().flatten()),
            **kwargs
        )

    # save a plot of the bandstrucuture
    if bandstructure_file is not None:
        ph.plot_bandstructure(phonon, Path(workdir) / bandstructure_file)


def initialize_phonopy_attach_calc(atoms, calc, supercell_matrix, displacement=0.01):
    """ phonopy preprocess returning supercells with attached calculator for FW """
    phonon, supercell, scs = ph.preprocess(atoms, supercell_matrix, displacement)
    for sc in scs:
        sc.calc = calc
    return phonon, supercell, scs
