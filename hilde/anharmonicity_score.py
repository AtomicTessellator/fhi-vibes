"""Functions to calculate anharmonicity scores as described in
    (future reference)"""
import numpy as np
import xarray as xr
from hilde.helpers.displacements import get_dR
from hilde.spglib.wrapper import get_symmetry_dataset


def _check_shape(f1, f2):
    """check if sizes of input data are compatible"""

    assert np.shape(f1) == np.shape(f2), (
        "Check shape of input arrays!: ",
        f1.shape,
        f2.shape,
    )


def get_r(in_f_data, in_f_model):
    r"""Calculate Pearson coefficient between f_data and f_model

    Refrence Website
    https://en.wikipedia.org/wiki/Pearson_correlation_coefficient

    Args:
        in_f_data: input data
        in_f_model: input model data

    Returns:
        float: Pearson coefficient
    """

    f_data = np.ravel(in_f_data)
    f_model = np.ravel(in_f_model)

    _check_shape(f_data, f_model)

    f_data_mean = np.mean(f_data, axis=0)
    f_model_mean = np.mean(f_model, axis=0)
    cov = (f_data - f_data_mean) @ (f_model - f_model_mean)
    sigma1 = np.sqrt((f_data - f_data_mean) @ (f_data - f_data_mean))
    sigma2 = np.sqrt((f_model - f_model_mean) @ (f_model - f_model_mean))

    return cov / sigma1 / sigma2


def get_r2(in_f_data, in_f_model):
    r"""Calculate coefficient of determination between f_data and f_model

    Refrence Website
    https://en.wikipedia.org/wiki/Coefficient_of_determination#Definitions

    Args:
        in_f_data: input data
        in_f_model: input model data

    Returns:
        float: Coefficient of Determination
    """

    f_data = np.ravel(in_f_data)
    f_model = np.ravel(in_f_model)

    _check_shape(f_data, f_model)

    f_data_mean = np.mean(f_data, axis=0)
    Sres = (f_data - f_model) @ (f_data - f_model)
    Stot = (f_data - f_data_mean) @ (f_data - f_data_mean)

    return 1 - Sres / Stot


def get_r2_per_sample(f_data, f_model, skip=0):
    """compute r2 for each sample"""
    data = []
    for fd, fm in zip(f_data[skip:], f_model[skip:]):
        data.append(get_r2(fd, fm))
    return np.array(data)


def get_r2_per_atom(
    forces_dft, forces_harmonic, ref_structure, reduce_by_symmetry=False
):
    """Compute r^2 score per atom in primitive cell. Optionally use symmetry.

    Args:
        forces_dft: forces from dft calculations
        forces_harmonic: forces from harmonic approximation
        ref_structure: reference structure for symmetry analysis
        reduce_by_symmetry: project on symmetry equivalent instead of primitive

    Returns:
        unique_atoms: the atoms from ref_structure for which r^2 was computed
        r2_per_atom: r^2 score for atoms in unique_atoms
    """

    sds = get_symmetry_dataset(ref_structure)

    if reduce_by_symmetry:
        compare_to = sds.equivalent_atoms
    else:
        compare_to = ref_structure.symbols

    if np.shape(forces_dft)[1] == 3 * len(ref_structure):
        compare_to = compare_to.repeat(3)

    unique_atoms, counts = np.unique(compare_to, return_counts=True)

    r2_atom = []
    for u in unique_atoms:
        # which atoms belong to the prototype?
        mask = compare_to == u
        # take the forces that belong to this atom
        f_dft = forces_dft[:, mask]
        f_ha = forces_harmonic[:, mask]
        # compute r^2
        r2_atom.append(get_r2(f_dft, f_ha))

    # tup = namedtuple("anharmonicity_score", "r2_per_atom unique counts")
    # ret = tup(np.array(r2_atom), unique_atoms, counts)

    attrs = {"unique_atoms": unique_atoms, "counts": counts}

    ret = xr.DataArray(np.array(r2_atom), attrs=attrs)

    return ret


def get_r2_MAE(in_f_data, in_f_model):
    r"""Calculate  1 - MAE/MA

    Args:
        in_f_data: input data
        in_f_model: input model data

    Returns:
        float: Coefficient of Determination
    """

    f_data = np.ravel(in_f_data)
    f_model = np.ravel(in_f_model)

    _check_shape(f_data, f_model)

    Sres = np.mean(abs(f_data - f_model))
    Stot = np.mean(abs(f_data))

    MAE = Sres / Stot

    return 1 - MAE


def get_forces_from_trajectory(trajectory, ref_structure=None, force_constants=None):
    """get forces from trajectory, consider to use `trajectory.set_forces_harmonic`

    Args:
        trajectory: list of Atoms objects with forces
        ref_structure: reference Atoms object
        force_constants: force constants in [3N, 3N] shape

    Returns:
        forces_dft: DFT forces in [N_steps, 3N] shape
        forces_harmonic: harmonic forces in [N_steps, 3N] shape
    """

    f_ha = get_forces_harmonic

    forces_dft = [a.get_forces().flatten() for a in trajectory]
    forces_ha = [f_ha(a, ref_structure, force_constants) for a in trajectory]

    return np.array(forces_dft), np.array(forces_ha)


def get_forces_harmonic(sc, ref_structure, force_constants):
    """helper function: compute forces from force_constants

    Args:
        sc: The distorted supercell
        ref_structure: The undistorted structure
        force_constants: The force constant matrix

    Returns:
        The harmonic forces
    """
    return -force_constants @ get_dR(sc, ref_structure).flatten()
