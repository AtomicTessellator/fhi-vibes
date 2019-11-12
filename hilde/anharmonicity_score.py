"""Functions to calculate anharmonicity scores as described in
    (future reference)"""
import numpy as np
import xarray as xr
from hilde.helpers import warn, progressbar
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


def get_r2(in_f_data, in_f_model, mean=True, silent=False):
    r"""Calculate coefficient of determination between f_data and f_model

    Refrence Website
    https://en.wikipedia.org/wiki/Coefficient_of_determination#Definitions

    Args:
        in_f_data: input data
        in_f_model: input model data
        mean (bool): take out mean of f_data
        silent (bool): make silent

    Returns:
        float: Coefficient of Determination
    """

    f_data = np.ravel(in_f_data)
    f_model = np.ravel(in_f_model)

    _check_shape(f_data, f_model)

    f_data_mean = np.mean(f_data, axis=0)

    # check f_data_mean, should be small
    if f_data_mean > 1e-1 and not silent:
        warn(f"|f_data.mean()| is {f_data_mean}", level=1)

    Sres = (f_data - f_model) @ (f_data - f_model)

    if mean:
        Stot = (f_data - f_data_mean) @ (f_data - f_data_mean)
    else:
        Stot = (f_data) @ (f_data)

    return 1 - Sres / Stot


def get_r2_per_sample(f_data, f_model, skip=0):
    """compute r2 for each sample"""
    data = []
    for fd, fm in zip(f_data[skip:], f_model[skip:]):
        data.append(get_r2(fd, fm))
    return np.array(data)


def get_r2_per_atom(
    forces_dft, forces_harmonic, ref_structure, reduce_by_symmetry=False, silent=True
):
    """Compute r^2 score per atom in primitive cell. Optionally use symmetry.

    Args:
        forces_dft: forces from dft calculations in shape [Nt, Na, 3]
        forces_harmonic: forces from harmonic approximation in shape [Nt, Na, 3]
        ref_structure: reference structure for symmetry analysis
        reduce_by_symmetry: project on symmetry equivalent instead of primitive
        silent (bool): silence the `get_r2` call

    Returns:
        unique_atoms: the atoms from ref_structure for which r^2 was computed
        r2_per_atom: r^2 score for atoms in unique_atoms
    """

    sds = get_symmetry_dataset(ref_structure)

    # check shape:
    if len(np.shape(forces_dft)) == 2:
        forces_dft = np.expand_dims(forces_dft, axis=0)
        forces_harmonic = np.expand_dims(forces_harmonic, axis=0)

    if reduce_by_symmetry:
        compare_to = sds.equivalent_atoms
    else:
        compare_to = ref_structure.numbers

    if np.shape(forces_dft)[1] == 3 * len(ref_structure):
        compare_to = compare_to.repeat(3)

    symbols = np.array(ref_structure.get_chemical_symbols())
    unique_atoms, counts = np.unique(compare_to, return_counts=True)

    r2_atom = []
    f_norms = []
    unique_symbols = []
    for u in unique_atoms:
        # which atoms belong to the prototype?
        mask = compare_to == u
        unique_symbols.append(symbols[mask][0])
        # take the forces that belong to this atom
        f_dft = forces_dft[:, mask]
        f_ha = forces_harmonic[:, mask]
        # compute r^2
        r2_atom.append(get_r2(f_dft, f_ha, silent=silent))
        f_norms.append(float(f_dft.std()))

    # compute weighted mean
    mean = (np.array(r2_atom) @ counts) / counts.sum()

    attrs = {
        "unique_atoms": unique_atoms,
        "symbols": unique_symbols,
        "counts": counts,
        "f_norms": f_norms,
        "mean": mean,
    }

    ret = xr.DataArray(np.array(r2_atom), attrs=attrs)

    return ret


def get_r2_MAE(in_f_data, in_f_model):
    r"""Calculate  1 - MAE/MA

    "r2 (RMA)": 1 - (1 - r2a) ** 2,

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

    r2a = 1 - MAE

    return 1 - (1 - r2a) ** 2


def get_r2_per_direction(f_data, f_model):
    """compute r2 for each Cartesian direction"""
    assert f_data.shape[-1] == 3
    directions = ("x", "y", "z")
    y = np.asarray(f_data).swapaxes(-1, 0)
    f = np.asarray(f_model).swapaxes(-1, 0)

    r2_direction = [get_r2(y[xx], f[xx]) for xx in range(len(directions))]
    return r2_direction


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


def get_dataframe(
    dataset, per_sample=False, per_direction=False, positive=False, by_symmetry=False
):
    """return anharmonicity dataframe for xarray.Dataset DS"""
    from pandas import DataFrame
    from hilde.helpers.converters import json2atoms

    DS = dataset

    ref_atoms = json2atoms(DS.attrs["reference atoms"])

    fs = DS.forces
    fhs = DS.forces_harmonic

    name = DS.attrs["System Name"]

    dct = {}
    if per_sample:
        for ii, (f, fh) in enumerate(zip(progressbar(fs), fhs)):

            key = f"{name}.{ii:04d}"
            dct[key] = {"r2 (RMS)": get_r2(f, fh)}

            r2s = get_r2_per_atom(f, fh, ref_atoms, reduce_by_symmetry=by_symmetry)
            dct[key].update({"r2_atom_mean": r2s.attrs["mean"]})
            d = {}
            for sym, r in zip(r2s.attrs["symbols"], r2s):  # enumerate(r2s):
                k = _key(sym, d)
                d[k] = float(r)

            dct[key].update(d)

            if per_direction:
                directions = ("x", "y", "z")
                r2s = get_r2_per_direction(f, fh)
                d = {}
                for ii, xx in enumerate(directions):
                    d[f"r2 [{xx}]"] = r2s[ii]

                dct[key].update(d)

    else:
        key = f"{name}"
        dct[key] = {"r2 (RMS)": get_r2(fs, fhs)}

        r2s = get_r2_per_atom(fs, fhs, ref_atoms, reduce_by_symmetry=by_symmetry)
        dct[key].update({"r2_atom_mean": r2s.attrs["mean"]})
        d = {}
        for sym, r in zip(r2s.attrs["symbols"], r2s):  # enumerate(r2s):
            k = _key(sym, d)
            d[k] = float(r)

        dct[key].update(d)

        if per_direction:
            directions = ("x", "y", "z")
            r2s = get_r2_per_direction(fs, fhs)
            d = {}
            for ii, xx in enumerate(directions):
                d[f"r2 [{xx}]"] = r2s[ii]

            dct[key].update(d)

    df = DataFrame(dct)
    # df[name] = df.mean(axis=1)
    return df.T


def _key(sym, d):
    key = f"r2 [{sym}]"
    if key in d:
        for ii in range(1, 100):
            key = f"r2 [{sym}{ii}]"
            if key not in d:
                break
    return key
