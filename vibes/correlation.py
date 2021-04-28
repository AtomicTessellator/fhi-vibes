import multiprocessing
from itertools import product

import numpy as np
import pandas as pd
import scipy.signal as sl
import xarray as xr
from vibes import dimensions, keys
from vibes.helpers import Timer, progressbar, talk
from vibes.helpers.warnings import warn


_prefix = "Correlation"
Timer.prefix = _prefix


def _talk(msg, verbose=True):
    return talk(msg, prefix=_prefix, verbose=verbose)


def _hann(nsamples: int):
    """Return one-side Hann function

    Args:
        nsamples (int): number of samples
    """
    return sl.windows.hann(2 * nsamples)[nsamples:]


hann = _hann


def get_unique_dimensions(dims: tuple) -> list:
    """make dimension labels unique by adding a counter to duplicate dims

    Example:
        ('time', 'a', 'a') -> ('time', 'a1', 'a2')

    """
    ar, index, counts = np.unique(dims, return_counts=True, return_inverse=True)
    ar = list(ar)

    # attach label to duplicate array elements
    for ii, count in enumerate(counts):
        if count > 1:
            elem = ar[ii]
            unique_elems = []
            for jj in range(count):
                unique_elem = elem + f"{jj+1}"
                unique_elems.append(unique_elem)

            ar[ii] = unique_elems

    # assemble new list and replace common unknowns
    new_dims = []
    for ii in index:
        elem = ar[ii]
        if isinstance(elem, list):
            elem = elem.pop(0)
            elem = dimensions.dim_replace.get(elem, elem)
        new_dims.append(elem)

    return new_dims


def _correlate(f1, f2, normalize=1, hann=True):
    """Compute correlation function for signal f1 and signal f2

    Reference:
        https://gitlab.com/flokno/python_recipes/-/blob/master/mathematics/
        correlation_function/autocorrelation.ipynb

    Args:
        f1: signal 1
        f2: signal 2
        normalize: no (0), by length (1), by lag (2)
        hann: apply Hann window
    Returns:
        the correlation function
    """
    a1, a2 = (np.asarray(f) for f in (f1, f2))
    Nt = min(len(a1), len(a2))

    if Nt != max(len(a1), len(a2)):
        msg = "The two signals are not equally long: "
        msg += f"len(a1), len(a2) = {len(a1)}, {len(a2)}"
        warn(msg)

    corr = sl.correlate(a1[:Nt], a2[:Nt])[Nt - 1 :]

    if normalize is True or normalize == 1:
        corr /= Nt
    elif normalize == 2:
        corr /= np.arange(Nt, 0, -1)

    if hann:
        corr *= _hann(Nt)

    return corr


# Frontends for _correlate
correlate = _correlate


def c1(X):
    """run `_correlate` with default parameters"""
    return correlate(X[0], X[1])


def get_autocorrelation(series, verbose=False, **kwargs):
    """Compute autocorrelation function of Series/DataArray

    Args:
        series ([N_t]): pandas.Series/xarray.DataArray with `time` axis in fs
        verbose (bool): be verbose
        kwargs: further arguments for `get_correlation`
    Return:
        DataArray ([N_t]): Autocorrelation function
    """
    timer = Timer("Compute autocorrelation function", verbose=verbose)

    autocorr = correlate(series, series, **kwargs)
    if isinstance(series, np.ndarray):
        result = autocorr
    elif isinstance(series, pd.Series):
        result = pd.Series(autocorr, index=series.index)
    elif isinstance(series, xr.DataArray):
        da = xr.DataArray(
            autocorr, dims=series.dims, coords=series.coords, name=keys.autocorrelation
        )
        result = da
    else:
        raise TypeError("`series` not of type ndarray, Series, or DataArray?")

    timer()

    return result


def _map_correlate(data, verbose=False):
    """compute autocorrelation with multiprocessing"""
    Nd = len(data)

    kw = {"len_it": Nd ** 2, "verbose": verbose}
    with multiprocessing.Pool() as p:
        res = [*progressbar(p.imap(c1, product(data, data), chunksize=Nd), **kw)]

    return np.asarray(res)


def get_autocorrelationNd(
    array: xr.DataArray,
    off_diagonal: bool = False,
    distribute: bool = True,
    verbose: bool = True,
    **kwargs,
) -> xr.DataArray:
    """compute velocity autocorrelation function for multi-dimensional xarray

    Default: Compute diagonal terms only (I, a, I, a)

    Args:
        array (xarray.DataArray [N_t, N_a, 3]): data
        off_diagonal_coords (bool): return off-diagonal coordinates (I, a, b)
        off_diagonal_atoms (bool): return off-diagonal atoms (I, J, a, b)
        kwargs: go to _correlate
    Returns:
        xarray.DataArray [N_t, N_a, 3]: autocorrelation along axis=0, or
        xarray.DataArray [N_t, N_a, 3, N_a, 3]: autocorrelation along axis=0
    """
    if not off_diagonal:
        return get_autocorrelation_diagonal(array)

    msg = "Get Nd autocorrelation"
    if off_diagonal:
        msg += " including off-diagonal terms"
    timer = Timer(msg, verbose=verbose)

    # memorize dimensions and shape of input array
    dims = array.dims
    Nt, *shape = np.shape(array)

    # compute autocorrelation including cross correlations
    # reshape to (X, Nt)
    data = np.asarray(array).reshape(Nt, -1).swapaxes(0, 1)

    # compute correlation function
    kw = {"verbose": verbose}
    if distribute:  # use multiprocessing
        corr = _map_correlate(data, **kw)
    else:
        corr = np.zeros([data.shape[0] ** 2, Nt])
        for ii, X in enumerate(progressbar([*product(data, data)], **kw)):
            corr[ii] = c1(X)

    # shape back and move time axis to the front
    corr = np.moveaxis(corr.reshape((*shape, *shape, Nt)), -1, 0)
    # create the respective dimensions
    new_dims = (dims[0], *dims[1:], *dims[1:])

    # make new dimensions unique
    new_dims = get_unique_dimensions(new_dims)
    # new_dims = np.unique(new_dims)

    df_corr = xr.DataArray(
        corr,
        dims=new_dims,
        coords=array.coords,
        name=f"{array.name}_{keys.autocorrelation}",
    )
    timer()

    return df_corr


def _exp(x, tau, y0):
    """compute exponential decay
        y = y0 * exp(-x / tau)
    """
    return y0 * np.exp(-x / tau)


def get_autocorrelation_diagonal(
    array: xr.DataArray, verbose: bool = True, **kwargs
) -> xr.DataArray:
    """compute autocorrelation function for each componetn in  multi-dimensional xarray

    Args:
        array [N_t, *dims]: data with time axis in the front
    Returns:
        xarray.DataArray [N_t, dims]: autocorrelation along axis=0

    """
    msg = f"Get autocorrelation for array of shape {array.shape}"
    timer = Timer(msg, verbose=verbose)

    # memorize dimensions and shape of input array
    Nt, *shape = np.shape(array)

    # compute autocorrelation for each dimension
    # move time axis to last index
    data = np.moveaxis(np.asarray(array), 0, -1)
    corr = np.zeros((*shape, Nt))
    for Ia in np.ndindex(*shape):
        tmp = correlate(data[Ia], data[Ia], **kwargs)
        corr[Ia] = tmp
    # move time axis back to front
    corr = np.moveaxis(corr, -1, 0)

    da_corr = array.copy()
    da_corr.data = corr
    da_corr.name = f"{array.name}_{keys.autocorrelation}"

    timer()

    return da_corr
