"""wrappers for scipy.integrate functions"""

import numpy as np
import pandas as pd
import xarray as xr
from scipy import integrate as si

from vibes import keys
from vibes.helpers import Timer, warn

_prefix = "Integration"
Timer.prefix = _prefix


def _cumtrapz(series, index=None, axis=0, initial=0):
    """wrap `scipy.integrate.cumtrapz`"""
    array = np.asarray(series)
    if index is not None and len(index) > 1:
        x = np.asarray(index)
    else:
        warn(f"index = {index}, use `x=None`", level=1)
        x = None
    ct = si.cumtrapz(array, x=x, axis=axis, initial=initial)

    return ct


def get_cumtrapz(series, **kwargs):
    """Compute cumulative trapezoid integral of ndarray, Series/DataArray

    Return:
        ndarray/Series/DataArray: cumulative trapezoid rule applied
    """
    if isinstance(series, np.ndarray):
        ctrapz = _cumtrapz(series, **kwargs)
        return ctrapz
    elif isinstance(series, pd.Series):
        ctrapz = _cumtrapz(series, index=series.index, **kwargs)
        return pd.Series(ctrapz, index=series.index)
    elif isinstance(series, xr.DataArray):
        ctrapz = _cumtrapz(series, index=series.coords, **kwargs)
        da = xr.DataArray(
            ctrapz, dims=series.dims, coords=series.coords, name=keys.cumtrapz
        )
        return da

    raise TypeError("`series` not of type ndarray, Series, or DataArray?")
