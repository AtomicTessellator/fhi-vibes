import multiprocessing
from itertools import product

import numpy as np
import pandas as pd
import scipy.optimize as so
import scipy.signal as sl
import xarray as xr

from vibes import keys
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


# Frontends for _correlate
correlate = _correlate


def c1(X):
    """run `_correlate` with default parameters"""
    return correlate(X[0], X[1])


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


def _map_correlate(data):
    """compute autocorrelation with multiprocessing"""
    Nd = len(data)

    with multiprocessing.Pool() as p:
        res = [
            *progressbar(p.imap(c1, product(data, data), chunksize=Nd), len_it=Nd ** 2,)
        ]

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
    msg = "Get Nd autocorrelation"
    if off_diagonal:
        msg += " including off-diagonal terms"
    timer = Timer(msg, verbose=verbose)

    # memorize dimensions and shape of input array
    dims = array.dims
    Nt, *shape = np.shape(array)

    # compute autocorrelation
    if off_diagonal:
        # reshape to (X, Nt)
        data = np.asarray(array).reshape(Nt, -1).swapaxes(0, 1)

        # compute correlation function
        if distribute:  # use multiprocessing
            corr = _map_correlate(data)
        else:
            corr = np.zeros([data.shape[0] ** 2, Nt])
            for ii, X in enumerate(progressbar([*product(data, data)])):
                corr[ii] = c1(X)

        # shape back and move time axis to the front
        corr = np.moveaxis(corr.reshape((*shape, *shape, Nt)), -1, 0)
        # create the respective dimensions
        new_dims = (dims[0], *dims[1:], *dims[1:])

    else:
        # move time axis to last index
        data = np.moveaxis(np.asarray(array), 0, -1)
        corr = np.zeros((*shape, Nt))
        for Ia in np.ndindex(*shape):
            tmp = correlate(data[Ia], data[Ia], **kwargs)
            corr[Ia] = tmp
        # move time axis back to front
        corr = np.moveaxis(corr, -1, 0)
        new_dims = dims

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


def get_correlation_time_estimate(
    series: pd.Series,
    tmin: float = 0.1,
    tmax: float = 5,
    tau0: float = 1,
    pre_smoothen_window: int = 100,
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
        pre_smoothen_window (int): use a window to pre-smoothen the data
        ps (bool): series.index given in ps (default: fs)
        return_series (bool): Return the corresponding exp. function as pd.Series
        verbose (bool): Be verbose.

    Returns:
        (float, float): tau (IN FS!!!), y0
    """
    kw = {"verbose": verbose, "prefix": _prefix}
    talk("Estimate correlation time from fitting exponential to normalized data", **kw)

    # smoothen the bare series
    series = series.rolling(pre_smoothen_window, min_periods=1).mean()
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
