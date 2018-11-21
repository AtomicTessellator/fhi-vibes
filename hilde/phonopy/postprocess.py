""" Provide a full highlevel phonopy workflow """

from pathlib import Path
import pickle
import numpy as np

from hilde.helpers.converters import dict2atoms
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import PhononRow
import hilde.phonopy.wrapper as ph
from hilde.phonopy import displacement_id_str
from hilde.structure.convert import to_Atoms
from hilde.trajectory import reader as traj_reader


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
