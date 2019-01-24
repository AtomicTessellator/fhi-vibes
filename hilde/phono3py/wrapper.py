"""
A leightweight wrapper for Phono3py
"""

from collections import namedtuple
import numpy as np
from phono3py.phonon3 import Phono3py
from hilde import konstanten as const
from hilde.helpers.supercell import make_cubic_supercell
from hilde.phonopy import enumerate_displacements
from hilde.structure.convert import to_Atoms, to_phonopy_atoms
from hilde.helpers.attribute_dict import AttributeDict as adict


defaults = adict(
    {
        "displacement": 0.03,
        "cutoff_pair_distance": 100.0,
        "symprec": 1e-5,
        "q_mesh": [11, 11, 11],
        "log_level": 2,
    }
)


def prepare_phono3py(
    atoms,
    supercell_matrix,
    fc3=None,
    phonon_supercell_matrix=None,
    fc2=None,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    q_mesh=defaults.q_mesh,
    displacement=defaults.displacement,
    symmetrize_fc3q=False,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
):
    """ Prepare a Phono3py object """

    ph_atoms = to_phonopy_atoms(atoms, wrap=True)

    supercell_matrix = get_3x3_matrix(supercell_matrix)

    if phonon_supercell_matrix is None:
        phonon_supercell_matrix = supercell_matrix
    else:
        phonon_supercell_matrix = get_3x3_matrix(phonon_supercell_matrix)

    phonon3 = Phono3py(
        ph_atoms,
        supercell_matrix=np.transpose(supercell_matrix),
        phonon_supercell_matrix=np.transpose(phonon_supercell_matrix),
        mesh=q_mesh,
        symprec=symprec,
        is_symmetry=True,
        symmetrize_fc3q=symmetrize_fc3q,
        frequency_factor_to_THz=const.omega_to_THz,
        log_level=log_level,
    )

    phonon3.generate_displacements(
        distance=displacement, cutoff_pair_distance=cutoff_pair_distance
    )

    if fc2 is not None:
        phonon3.set_fc2(fc2)
    if fc3 is not None:
        phonon3.set_fc3(fc3)

    return phonon3


def preprocess(
    atoms,
    supercell_matrix,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    q_mesh=defaults.q_mesh,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
):
    """
    Set up a Phono3py object and generate all the supercells necessary for the 3rd order
    """

    phonon3 = prepare_phono3py(
        atoms,
        supercell_matrix=supercell_matrix,
        cutoff_pair_distance=cutoff_pair_distance,
        q_mesh=q_mesh,
        displacement=displacement,
        symprec=symprec,
        log_level=log_level,
    )

    supercell = to_Atoms(
        phonon3.get_supercell(),
        info={
            "supercell": True,
            "supercell_matrix": phonon3.get_supercell_matrix().T.flatten().tolist(),
        },
    )

    scells = phonon3.get_supercells_with_displacements()
    supercells_with_disps = [to_Atoms(cell) for cell in scells]

    enumerate_displacements(supercells_with_disps)

    pp = namedtuple(
        "phono3py_preprocess", "phonon supercell supercells_with_displacements"
    )

    return pp(phonon3, supercell, supercells_with_disps)
