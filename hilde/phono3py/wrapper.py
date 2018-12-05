"""
A leightweight wrapper for Phono3py
"""

from collections import namedtuple
import numpy as np
from phono3py.phonon3 import Phono3py
from hilde import konstanten as const
from hilde.phonopy import enumerate_displacements
from hilde.structure.convert import to_Atoms, to_phonopy_atoms
from hilde.helpers.config import AttributeDict as adict


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
    fc2_supercell_matrix=[1,1,1],
    fc3_supercell_matrix=[1,1,1],
    q_mesh=defaults.q_mesh,
    fc2=None,
    fc3=None,
    disp=defaults.displacement,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    symmetrize_fc3q=False,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
):
    """ Prepare a Phono3py object """

    ph_atoms = to_phonopy_atoms(atoms, wrap=True)

    phonon3 = Phono3py(
        ph_atoms,
        supercell_matrix=fc3_supercell_matrix,
        phonon_supercell_matrix=fc2_supercell_matrix,
        mesh=q_mesh,
        symprec=symprec,
        is_symmetry=True,
        symmetrize_fc3q=symmetrize_fc3q,
        frequency_factor_to_THz=const.eV_to_THz,
        log_level=log_level,
    )

    phonon3.generate_displacements(
        distance=disp, cutoff_pair_distance=cutoff_pair_distance
    )

    if fc2 is not None:
        phonon3.set_fc2(fc2)
    if fc3 is not None:
        phonon3.set_fc3(fc3)

    return phonon3


def preprocess(
    atoms,
    fc2_supercell_matrix,
    fc3_supercell_matrix,
    q_mesh=defaults.q_mesh,
    disp=defaults.displacement,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
):
    """ Set up a Phono3py object and generate all the necessary supercells """

    phonon3 = prepare_phono3py(
        atoms,
        fc2_supercell_matrix=np.array(fc2_supercell_matrix).reshape(3, 3),
        fc3_supercell_matrix=np.array(fc3_supercell_matrix).reshape(3, 3),
        q_mesh=q_mesh,
        disp=disp,
        cutoff_pair_distance=cutoff_pair_distance,
        symprec=symprec,
        log_level=log_level,
    )

    # phonpoy supercells
    fc2_supercell = to_Atoms(
        phonon3.get_phonon_supercell(),
        info={
            "supercell": True,
            "supercell_matrix": np.asarray(fc2_supercell_matrix).flatten().tolist(),
        },
    )

    fc3_supercell = to_Atoms(
        phonon3.get_supercell(),
        info={
            "supercell": True,
            "supercell_matrix": np.asarray(fc3_supercell_matrix).flatten().tolist(),
        },
    )

    scells = phonon3.get_phonon_supercells_with_displacements()
    fc2_supercells_with_disps = [to_Atoms(cell) for cell in scells]

    scells = phonon3.get_supercells_with_displacements()
    fc3_supercells_with_disps = [to_Atoms(cell) for cell in scells]

    enumerate_displacements([*fc2_supercells_with_disps, *fc3_supercells_with_disps])

    cols = (
        "phonon3",
        "fc2_supercell",
        "fc3_supercell",
        "fc2_supercells_with_displacements",
        "fc3_supercells_with_displacements",
    )
    pp = namedtuple("phono3py_preprocess", " ".join(cols))

    return pp(
        phonon3,
        fc2_supercell,
        fc3_supercell,
        fc2_supercells_with_disps,
        fc3_supercells_with_disps,
    )


def get_forces(supercells_computed):
    """ Return force_sets taking care of supercells that were not computed
    because of cutoff. """

    zero_force = np.zeros([len(supercells_computed[0]), 3])
    force_sets = []
    for scell in supercells_computed:
        if scell is None:
            force_sets.append(zero_force)
        else:
            force_sets.append(scell.get_forces())

    if len(force_sets) != len(supercells_computed):
        print(
            "len(force_sets), len(supercells):",
            len(force_sets),
            len(supercells_computed),
        )
        raise RuntimeError("Number of computed supercells incorrect.")

    return force_sets
