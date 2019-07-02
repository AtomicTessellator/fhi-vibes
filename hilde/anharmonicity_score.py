"""Functions to calculate anharmonicity scores as described in
    (future reference)"""

import numpy as np
from hilde.harmonic_analysis.displacements import get_dR
from hilde.spglib.wrapper import get_symmetry_dataset


def get_r2(in_f_data, in_f_model):
    r"""Calculate coefficient of determination between f_data and f_model

    Refrence Website
    https://en.wikipedia.org/wiki/Coefficient_of_determination#Definitions

    Parameters
    -----------
    in_f_data: array
        input data
    in_f_model: array
        input model data

    Returns
    -------
    r2: float
        Coefficient of Determination
    """

    f_data = np.ravel(in_f_data)
    f_model = np.ravel(in_f_model)

    f_data_mean = np.mean(f_data, axis=0)
    Sres = (f_data - f_model) @ (f_data - f_model)
    Stot = (f_data - f_data_mean) @ (f_data - f_data_mean)

    return 1 - Sres / Stot


def get_r2_per_atom(
    forces_dft, forces_harmonic, ref_structure, reduce_by_symmetry=False
):
    """Compute r^2 score per atom in primitive cell. Optionally use symmetry.

    Parameters
    ----------
    forces_dft: list
        forces from dft calculations
    forces_harmonic: list
        forces from harmonic approximation
    ref_structure: ase.atoms.Atoms
        reference structure for symmetry analysis
    reduce_by_symmetry: bool
        project on symmetry equivalent instead of primitive

    Returns
    -------
    unique_atoms: list
        the atoms from ref_structure for which r^2 was computed
    r2_per_atom: list
        r^2 score for atoms in unique_atoms
    """

    sds = get_symmetry_dataset(ref_structure)

    if reduce_by_symmetry:
        compare_to = sds.equivalent_atoms
    else:
        compare_to = sds.mapping_to_primitive

    unique_atoms = np.unique(compare_to)

    r2_atom = []
    for u in unique_atoms:
        # which atoms belong to the prototype?
        mask = compare_to.repeat(3) == u
        # take the forces that belong to this atom
        f_dft = forces_dft[:, mask]
        f_ha = forces_harmonic[:, mask]
        # compute r^2
        r2_atom.append(get_r2(f_dft, f_ha))

    return unique_atoms, r2_atom


def get_forces_from_trajectory(trajectory, ref_structure=None, force_constants=None):
    """get forces from trajectory

    Parameters
    ----------
    trajectory: list
        list of Atoms objects with forces
    ref_structure: ase.atoms.Atoms
        reference Atoms object
    force_constants: np.ndarray
        force constants in [3N, 3N] shape

    Returns
    -------
    forces_dft: np.ndarray
        DFT forces in [N_steps, 3N] shape
    forces_harmonic: np.ndarray
        harmonic forces in [N_steps, 3N] shape
    """

    f_ha = get_harmonic_forces

    forces_dft = [a.get_forces().flatten() for a in trajectory]
    forces_ha = [f_ha(a, ref_structure, force_constants) for a in trajectory]

    return np.array(forces_dft), np.array(forces_ha)


def get_harmonic_forces(sc, ref_structure, force_constants):
    """helper function: compute forces from force_constants

    Parameters
    ----------
    sc: ase.atoms.Atoms
        The distorted supercell
    ref_structure: ase.atoms.Atoms
        The undistorted structure
    force_constants: np.ndarray
        The force constant matrix

    Results
    -------
    np.ndarray
        The harmonic forces
    """
    return -force_constants @ get_dR(sc, ref_structure).flatten()
