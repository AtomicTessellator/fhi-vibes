"""
A leightweight wrapper for Phono3py
"""

import numpy as np
from phono3py.phonon3 import Phono3py
from hilde import konstanten as const
from hilde.phonopy import get_supercells_with_displacements
from hilde.structure.convert import to_phonopy_atoms
from hilde.helpers.attribute_dict import AttributeDict as adict
from hilde.helpers.numerics import get_3x3_matrix

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
    displacement_dataset=None,
    q_mesh=defaults.q_mesh,
    displacement=defaults.displacement,
    symmetrize_fc3q=False,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
    **kwargs,
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

    if displacement_dataset:
        phonon3.set_displacement_dataset(displacement_dataset)

    phonon3.generate_displacements(
        distance=displacement, cutoff_pair_distance=cutoff_pair_distance
    )

    if fc2 is not None:
        phonon3.set_fc2(fc2)
    if fc3 is not None:
        phonon3.set_fc3(fc3)

    return phonon3

def preprocess_fireworks(
    atoms,
    calc,
    kpt_density=None,
    supercell_matrix=None,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    up_kpoint_from_pc=False,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
):
    """ wrap the Phono3py preprocess for workflow """
    if supercell_matrix is None:
        raise InputError("The supercell_matrix must be defined")

    # Phonopy preprocess
    phonon3, _, supercell, _, scs = phono3py_preprocess(
        atoms=atoms,
        fc3_supercell_matrix=supercell_matrix,
        n_atoms_in_sc_3=n_atoms_in_sc,
        disp=displacement,
        cutoff_pair_distance=cutoff_pair_distance,
        symprec=symprec,
    )
    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    if kpt_density is not None:
        update_k_grid(supercell, calc, kpt_density)

    metadata = metadata2dict(atoms, calc, phonon3)

    to_run_scs = []
    for sc in scs:
        if sc:
            sc.calc = calc
            to_run_scs.append(sc)
    return calc, supercell, to_run_scs, phonon3, metadata

def preprocess(
    atoms,
    supercell_matrix,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    q_mesh=defaults.q_mesh,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
    **kwargs,
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

    return get_supercells_with_displacements(phonon3)
