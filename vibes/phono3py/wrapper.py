""" A leightweight wrapper for Phono3py """

import numpy as np
from phono3py import Phono3py

from vibes import konstanten as const
from vibes.helpers.numerics import get_3x3_matrix
from vibes.phonopy import get_supercells_with_displacements
from vibes.structure.convert import to_phonopy_atoms

from ._defaults import defaults


def prepare_phono3py(
    atoms,
    supercell_matrix,
    fc2=None,
    fc3=None,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    displacement_dataset=None,
    is_diagonal=defaults.is_diagonal,
    q_mesh=defaults.q_mesh,
    displacement=defaults.displacement,
    symmetrize_fc3q=False,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
    **kwargs,
):
    """Prepare a Phono3py object

    Args:
        atoms: ase.atoms.Atoms
        supercell_matrix: np.ndarray
        fc2: np.ndarray
        fc3: np.ndarray
        cutoff_pair_distance: float
        displacement_dataset: dict
        is_diagonal: bool
        mesh: np.ndarray
        displacement: float
        symmetrize_fc3q: bool
        symprec: float
        log_level: int

    Returns:
        phono3py.Phono3py
    """

    ph_atoms = to_phonopy_atoms(atoms, wrap=True)

    supercell_matrix = get_3x3_matrix(supercell_matrix)

    phonon3 = Phono3py(
        ph_atoms,
        supercell_matrix=np.transpose(supercell_matrix),
        mesh=q_mesh,
        symprec=symprec,
        is_symmetry=True,
        symmetrize_fc3q=symmetrize_fc3q,
        frequency_factor_to_THz=const.omega_to_THz,
        log_level=log_level,
    )

    if displacement_dataset is not None:
        phonon3.set_displacement_dataset(displacement_dataset)

    phonon3.generate_displacements(
        distance=displacement,
        cutoff_pair_distance=cutoff_pair_distance,
        is_diagonal=is_diagonal,
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
    is_diagonal=defaults.is_diagonal,
    q_mesh=defaults.q_mesh,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
    **kwargs,
):
    """Set up a Phono3py object and generate all the supercells necessary for the 3rd order

    Args:
        atoms: ase.atoms.Atoms
        supercell_matrix: np.ndarray
        cutoff_pair_distance: float
        is_diagonal: bool
        q_mesh: np.ndarray
        displacement: float
        symprec: float
        log_level: int

    Returns:
        phonon3: phono3py.Phono3py
        supercell: ase.atoms.Atoms
        supercells_with_disps: list of ase.atoms.Atoms
    """

    phonon3 = prepare_phono3py(
        atoms,
        supercell_matrix=supercell_matrix,
        cutoff_pair_distance=cutoff_pair_distance,
        is_diagonal=is_diagonal,
        q_mesh=q_mesh,
        displacement=displacement,
        symprec=symprec,
        log_level=log_level,
    )

    return get_supercells_with_displacements(phonon3)
