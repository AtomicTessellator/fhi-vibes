""" A leightweight wrapper for Phono3py """

import numpy as np
from phono3py.phonon3 import Phono3py
from hilde import konstanten as const
from hilde.phonopy import get_supercells_with_displacements
from hilde.structure.convert import to_phonopy_atoms
from hilde.helpers.numerics import get_3x3_matrix
from ._defaults import defaults


def prepare_phono3py(
    atoms,
    supercell_matrix,
    fc3=None,
    phonon_supercell_matrix=None,
    fc2=None,
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

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        primitive cell for the calculation
    supercell_matrix: np.ndarray
        supercell matrix for the third order phonons
    fc3: np.ndarray
        Third order force constant matrix
    phonon_supercell_matrix: np.ndarray
        supercell matrix for the second order phonons
    fc2: np.ndarray
        second order force constant matrix
    cutoff_pair_distance: float
        All pairs further apart than this cutoff are ignored
    displacement_dataset: dict
        The displacement_dataset for the third order phonons
    is_diagonal: bool
        Whether allow diagonal displacements of Atom 2 or not
    q_mesh: np.ndarray
        q-point interpolation mesh postprocessing
    displacement: float
        magnitude of the displacement
    symmetrize_fc3q: bool
        If True symmetrize the third order interactions
    symprec: float
        distance tolerance for determining the sapce group/symmetry
    log_level: int
        How much information should be streamed to the console

    Returns
    -------
    phonon3: phono3py.phonon3.Phono3py
        The Phono3py object for the calculation
    """

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

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        primitive cell for the calculation
    supercell_matrix: np.ndarray
        supercell matrix for the third order phonons
    cutoff_pair_distance: float
        All pairs further apart than this cutoff are ignored
    is_diagonal: bool
        Whether allow diagonal displacements of Atom 2 or not
    q_mesh: np.ndarray
        q-point interpolation mesh postprocessing
    displacement: float
        magnitude of the displacement
    symprec: float
        distance tolerance for determining the sapce group/symmetry
    log_level: int
        How much information should be streamed to the console

    Returns
    -------
    phonon3: phono3py.phonon3.Phono3py
        The Phono3py object with displacement_dataset, and displaced supercells
    supercell: ase.atoms.Atoms
        The undisplaced supercell
    supercells_with_disps: list of ase.atoms.Atoms
        All of the supercells with displacements
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
