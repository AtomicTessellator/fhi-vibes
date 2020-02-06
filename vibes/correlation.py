import numpy as np
import pandas as pd
import scipy.signal as sl
import xarray as xr

from vibes import keys
from vibes.helpers import Timer
from vibes.helpers.warnings import warn

_prefix = "Correlation"
Timer.prefix = _prefix


def _correlate(f1, f2, normalize=True, window=True):
    """Compute correlation function for signal f1 and signal f2

    Reference:
        https://gitlab.com/flokno/python_recipes/-/blob/master/mathematics/
        correlation_function/autocorrelation.ipynb

    Args:
        f1: signal 1
        f2: signal 2
        normalize: normalize correlation function
        window: apply Hann window
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

    if normalize:
        corr /= np.arange(Nt, 0, -1)

    if window:
        corr *= sl.windows.hann(2 * Nt)[Nt:]

    return corr


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

    autocorr = _correlate(series, series, **kwargs)
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
