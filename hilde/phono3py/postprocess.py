""" Provide a full highlevel phonopy workflow """

from pathlib import Path
import pickle
import numpy as np

from phono3py.file_IO import write_fc3_to_hdf5
from hilde.helpers.converters import dict2atoms
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import PhononRow
from hilde.phonopy import displacement_id_str
from hilde.trajectory import reader as traj_reader


def postprocess(
    phonon3,
    symmetrize_fc3=False,
    workdir=".",
    trajectory="trajectory.yaml",
    force_constants_file="force_constants_3.hdf5",
    pickle_file="phono3py.pick",
    **kwargs,
):
    """ Phonopy postprocess """

    trajectory = Path(workdir) / trajectory

    if trajectory.is_file():
        calculated_atoms = traj_reader(trajectory)
    else:
        raise ValueError("Either calculated_atoms or trajectory must be defined")

    forces_shape = calculated_atoms[0].get_forces().shape

    force_sets_fc3 = []
    used_forces = 0
    for ii, cell in phonon3.get_supercells_with_displacements():
        if cell is not None:
            ref_atoms = calculated_atoms[used_forces]
            if ref_atoms.info[displacement_id_str] == ii:
                force_sets_fc3.append(ref_atoms.get_forces())
                used_forces += 1
        else:
            force_sets_fc3.append(np.zeros(forces_shape))

    # compute force constants
    if symmetrize_fc3:
        phonon3.produce_fc3(
            force_sets_fc3,
            is_translational_symmetry=True,
            is_permutation_symmetry=True,
            is_permutation_symmetry_fc2=True,
        )

    else:
        phonon3.produce_fc3(force_sets_fc3)

    write_fc3_to_hdf5(phonon3.get_fc3(), filename=force_constants_file)

    exit(
        f"*** Force Sets provided in {pickle_file}, "
        "full postprocess not yet supported"
    )
