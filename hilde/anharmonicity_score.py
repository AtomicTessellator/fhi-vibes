"""Functions to calculate anharmonicity scores as described in
    (future reference)"""

import numpy as np
from hilde.trajectory import reader
from hilde.harmonic_analysis.displacements import get_dR
from hilde.helpers.converters import dict2atoms
from hilde.spglib.wrapper import get_symmetry_dataset


def get_r2(in_f_data, in_f_model):
    r"""Calculate coefficient of determination between f_data and f_model

    Reference: https://en.wikipedia.org/wiki/Coefficient_of_determination#Definitions

    Args:
        in_f_data (array): input data
        in_f_model (array): input model data

    Returns:
        r2 (float): Coefficient of Determination

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

    Args:
        forces_dft (list): forces from dft calculations
        forces_harmonic (list): forces from harmonic approximation
        ref_structure (ase.Atoms): reference structure for symmetry analysis
        reduce_by_symmetry (bool): project on symmetry equivalent instead of primitive

    Returns:
        unique_atoms (list): the atoms from ref_structure for which r^2 was computed
        r2_per_atom (list): r^2 score for atoms in unique_atoms

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

    Args:
        trajectory (list): list of ase.Atoms objects with forces
        ref_structure (ase.Atoms): reference Atoms object
        force_constants (np.ndarray): force constants in [3N, 3N] shape

    Returns:
        forces_dft (np.ndarray): DFT forces in [N_steps, 3N] shape
        forces_harmonic (np.ndarray): harmonic forces in [N_steps, 3N] shape
    """

    f_ha = get_harmonic_forces

    forces_dft = [a.get_forces().flatten() for a in trajectory]
    forces_ha = [f_ha(a, ref_structure, force_constants) for a in trajectory]

    return np.array(forces_dft), np.array(forces_ha)


def get_harmonic_forces(sc, ref_structure, force_constants):
    """helper function: compute forces from force_constants"""
    return -force_constants @ get_dR(sc, ref_structure).flatten()


def get_force_sets_from_trajectory(
    trajectory="statistical_sampling/trajectory.son",
    force_constants=None,
    ref_structure=None,
    temperature=None,
):
    # Get the displaced atoms
    scs_displaced, meta = reader(trajectory, True)

    if ref_structure is None:
        ref_structure = dict2atoms(meta["atoms"])

    if force_constants is None:
        force_constants = np.array(meta["force_constants"])

    forces_dft = dict()
    forces_harmonic = dict()
    for sc in scs_displaced:
        if temperature and sc.info["temperature"] != temperature:
            continue
        if sc.info["temperature"] not in forces_dft:
            forces_dft[sc.info["temperature"]] = list()
            forces_harmonic[sc.info["temperature"]] = list()
        forces_dft[sc.info["temperature"]] += list(sc.get_forces().flatten())
        forces_harmonic[sc.info["temperature"]] += list(
            get_harmonic_forces(sc, ref_structure, force_constants)
        )
    if temperature is None:
        return forces_dft, forces_harmonic
    return forces_dft[temperature], forces_harmonic[temperature]


def reshape_forces(forces_dft, forces_harmonic):
    forces_dft = forces_dft.copy().reshape(-1, 3)
    forces_harmonic = forces_harmonic.copy().reshape(-1, 3)
    return forces_dft, forces_harmonic


def get_deviation(forces_dft, forces_harmonic, ref_structure):
    forces_dft, forces_harmonic = reshape_forces(forces_dft, forces_harmonic)

    vol = ref_structure.get_volume()

    diff = np.abs(forces_dft - forces_harmonic)
    dev_vol_nromalized = np.linalg.norm(diff / vol, axis=1)

    for ii in range(3):
        diff[:, ii] /= np.linalg.norm(forces_harmonic, axis=1)

    dev_force_normalized = diff.copy()

    return dev_vol_nromalized, dev_force_normalized


def get_rmse(forces_dft, forces_harmonic):
    rmse = np.sqrt(
        np.sum((forces_dft - forces_harmonic) ** 2.0) / forces_harmonic.shape[0]
    )
    return rmse


def get_normalized_rmse(forces_dft, forces_harmonic, ref_structure):
    vol = ref_structure.get_volume()
    rmse_vol_normalized = get_rmse(forces_dft, forces_harmonic) / vol

    forces_dft, forces_harmonic = reshape_forces(forces_dft, forces_harmonic)
    force_harmonic_mag = np.linalg.norm(forces_harmonic, axis=1)
    forces_dft /= force_harmonic_mag
    forces_harmonic /= force_harmonic_mag
    rmse_force_normalized = get_rmse(forces_dft.flatten(), forces_harmonic.flatten())

    return rmse_vol_normalized, rmse_force_normalized


def get_anharmonicity_score(
    trajectory="statistical_sampling/trajectory.son",
    force_constants=None,
    ref_structure=None,
    r2_file="r2.dat",
    **kwargs
):
    """
    Calculates the R^2 anharmonicity parameter for a given structure
    Args:
        trajectory (str): File name for the trajectory of statistically sampled supercell displacements
        phonon_file (str): File name for the phonopy trajectory used to generate the force constants used to generate the structures in trajectory
        r2_file(str): File name to output the temperature and r^2 value
    Returns (tuple(float, float)): Temperature, Coefficient of Determination between actual forces and harmonic ones
    """
    r2 = []
    forces_dft, forces_harmonic = get_force_sets_from_trajectory(
        trajectory, force_constants, ref_structure
    )
    for temp in sorted(forces_dft.keys()):
        r2.append(
            [temp, get_r2(np.array(forces_dft[temp]), np.array(forces_harmonic[temp]))]
        )
    np.savetxt(r2_file, r2, header="Temperature, R^2")
    return r2
