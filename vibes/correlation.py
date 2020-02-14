from itertools import chain
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.optimize as so
import scipy.signal as sl
import xarray as xr

from vibes import dimensions, keys
from vibes.helpers import Timer, progressbar, talk
from vibes.helpers.warnings import warn

_prefix = "Correlation"
Timer.prefix = _prefix


def _hann(nsamples: int):
    """Return one-side Hann function

    Args:
        nsamples (int): number of samples
    """
    return sl.windows.hann(2 * nsamples)[nsamples:]


hann = _hann


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


correlate = _correlate


def get_autocorrelation(series, verbose=True, **kwargs):
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


def get_autocorrelationNd(
    array: xr.DataArray,
    off_diagonal_coords: bool = False,
    off_diagonal_atoms: bool = False,
    cache: bool = True,
    verbose: bool = True,
    **kwargs,
) -> xr.DataArray:
    """compute velocity autocorrelation function for multi-dimensional xarray

    Default: Compute diagonal terms only (I, a, I, a)

    Args:
        array (xarray.DataArray [N_t, N_a, 3]): data
        off_diagonal_coords (bool): return off-diagonal coordinates (I, a, b)
        off_diagonal_atoms (bool): return off-diagonal atoms (I, J, a, b)
        cache: cache
        kwargs: go to _correlate
    Returns:
        xarray.DataArray [N_t, N_a, 3]: autocorrelation along axis=0, or
        xarray.DataArray [N_t, N_a, 3, 3]: autocorrelation along axis=0, or
        xarray.DataArray [N_t, N_a, N_a, 3, 3]: autocorrelation along axis=0
    """
    msg = "Get Nd autocorrelation"
    if off_diagonal_atoms:
        off_diagonal_coords = True
        msg += " including off-diagonal terms and coordinates"
    timer = Timer(msg, verbose=verbose)

    Nt, *shape = np.shape(array)

    dims = array.dims

    if dims == dimensions.time_vec:
        new_dims = dimensions.time_tensor
    elif dims == dimensions.time_atom_vec:
        new_dims = dimensions.time_atom_tensor
    else:
        raise TypeError(f"FIXME: Dimensions not implemented: {dims}")

    # move time axis to end
    data = np.moveaxis(np.asarray(array), 0, -1)

    # expand fake atoms axis?
    if len(data.shape) < 3:
        if off_diagonal_atoms:
            raise ValueError(f"Presumably not atoms given, inspect! Array dims:", dims)
        shape = (1, *shape)
        data = data[np.newaxis, :]

    # compute autocorrelation
    if off_diagonal_atoms and off_diagonal_coords:
        new_dims = dimensions.time_atom_atom_tensor
        corr = np.zeros([*np.repeat(shape, 2), Nt])
        for Ia in progressbar([*np.ndindex(*shape)]):
            for Jb in np.ndindex(*shape):
                tmp = correlate(data[Ia], data[Jb], **kwargs)
                idx = tuple(chain.from_iterable(zip(Ia, Jb)))  # flatten
                corr[idx] = tmp
    elif off_diagonal_coords:
        corr = np.zeros([*shape, shape[-1], Nt])
        for Ia in np.ndindex(*shape):
            for b in range(shape[-1]):
                Jb = (Ia[0], b)
                tmp = correlate(data[Ia], data[Jb], **kwargs)
                corr[(*Ia, b)] = tmp
    else:
        new_dims = new_dims[:-1]
        corr = np.zeros([*shape, Nt])
        for Ia in np.ndindex(*shape):
            tmp = correlate(data[Ia], data[Ia], **kwargs)
            corr[Ia] = tmp

    # move time axis back to 0 and remove aux. axis
    corr = np.moveaxis(corr, -1, 0).squeeze()

    df_corr = xr.DataArray(
        corr,
        dims=new_dims,
        coords=array.coords,
        name=f"{array.name}_{keys.autocorrelation}",
    )

    if cache:
        outfile = Path(f"{keys.cache}") / f"{get_autocorrelationNd.__name__}.nc"
        outfile.parent.mkdir(exist_ok=True)
        df_corr.to_netcdf(outfile)
        talk(f".. cache to {outfile}")

    timer()

    return df_corr


def _exp(x, tau, y0):
    """compute exponential decay
        y = y0 * exp(-x / tau)
    """
    return y0 * np.exp(-x / tau)


def get_correlation_time_estimate(
    series: pd.Series,
    tmin: float = 0.1,
    tmax: float = 5,
    tau0: float = 1,
    # pre_smoothen_window: int = 10,
    ps: bool = False,
    return_series: bool = False,
    verbose: bool = True,
) -> (float, float):
    """estimate correlation time of series by fitting an exponential its head to

        y = y0 * exp(-x / tau)

    Args:
        series (pd.Series): the time series
        tmin (float, optional): Start fit (ps). Defaults to 0.1.
        tmax (float, optional): End fit (ps). Defaults to 5.
        tau0 (float, optional): Guess for correlation time. Defaults to 1.
        # pre_smoothen_window (int): use a window to pre-smoothen the data
        ps (bool): series.index given in ps (default: fs)
        return_series (bool): Return the corresponding exp. function as pd.Series
        verbose (bool): Be verbose.

    Returns:
        (float, float): tau (IN FS!!!), y0
    """
    kw = {"verbose": verbose, "prefix": _prefix}
    talk("Estimate correlation time from fitting exponential to normalized data", **kw)

    # smoothen the bare series
    # series = series.rolling(pre_smoothen_window, min_periods=0).mean()
    series = series.copy().dropna()

    if not ps:
        series.index /= 1000

    y = series[tmin:tmax] / series.iloc[0]

    if len(y) < 10:
        warn("Not enough data to produce estimate, return default values")
        return tau0, 1

    # to ps
    x = y.index
    bounds = ([0, 0], [20 * tau0, 0.5])

    (tau, y0), _ = so.curve_fit(_exp, x, y, bounds=bounds)

    talk(f"Correlation time tau:       {tau:7.3f} ps", **kw)
    talk(f"Intersection y0:            {y0:7.3f}", **kw)
    talk(f"90% integral (2.30 * tau):  {2.3 * tau:7.3f} ps", **kw)
    talk(f"95% integral (3.00 * tau):  {3.0 * tau:7.3f} ps", **kw)
    talk(f"99% integral (4.61 * tau):  {4.605 * tau:7.3f} ps", **kw)

    if not ps:
        tau *= 1000

    return tau, y0


def get_autocorrelation_exponential(
    series: pd.Series, ps: bool = False, verbose: bool = True, **kwargs
) -> pd.Series:
    """Return exponential fit to time series assuming an exp. decay

    Args:
        series (pd.Series): the time series (time in fs)
        ps (bool): series.index given in ps instead of fs
        verbose (bool, optional): Be verbose. Defaults to True.
        kwargs (dict): kwargs for get_correlation_time_estimate

    Returns:
        pd.Series
    """
    tau, y0 = get_correlation_time_estimate(series, ps=ps, verbose=verbose)

    if ps:
        return pd.Series(_exp(series.index, tau, y0), index=series.index)
    else:
        return pd.Series(_exp(series.index, tau * 1000, y0), index=series.index)


# @numba.jit(parallel=True)
# def _autocorrelationNd(
#     data: np.ndarray, corr: np.ndarray, Nt: int, shape: list
# ) -> np.ndarray:
#     for Ia, d1 in np.ndenumerate(data):
#         for Jb, d2 in np.ndenumerate(data):
#             tmp = sl.correlate(d1, d2)[Nt - 1 :]
#             # tmp = _correlate(d1, d2)
#             idx = tuple(chain.from_iterable(zip(Ia, Jb)))  # flatten
#             corr[idx] = tmp
#
#     return corr
