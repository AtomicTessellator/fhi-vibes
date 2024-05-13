"""wrappers for scipy.integrate functions"""

import numpy as np
import pandas as pd
import xarray as xr
from scipy import integrate as si

from vibes import keys
from vibes.helpers import Timer, warn

_prefix = "Integration"
Timer.prefix = _prefix


def trapz(series, index=None, axis=0, initial=0):
    """Wrap `scipy.integrate.trapz`"""
    array = np.asarray(series)
    if index is not None and len(index) > 1:
        x = np.asarray(index)
    else:
        warn(f"index = {index}, use `x=None`", level=1)
        x = None
    return si.trapz(array, x=x, axis=axis)


def cumtrapz(series, index=None, axis=0, initial=0):
    """Wrap `scipy.integrate.cumtrapz`"""
    array = np.asarray(series)
    if index is not None and len(index) > 1:
        x = np.asarray(index)
    else:
        warn(f"index = {index}, use `x=None`", level=1)
        x = None
    return si.cumtrapz(array, x=x, axis=axis, initial=initial)


def get_cumtrapz(series, **kwargs):
    """
    Compute cumulative trapezoid integral of ndarray, Series/DataArray

    Return:
    ------
        ndarray/Series/DataArray: cumulative trapezoid rule applied

    """
    if isinstance(series, np.ndarray):
        return cumtrapz(series, **kwargs)

    if isinstance(series, pd.Series):
        ctrapz = cumtrapz(series, index=series.index, **kwargs)
        return pd.Series(ctrapz, index=series.index)

    if isinstance(series, xr.DataArray):
        try:
            index = np.asarray(series[keys.time])
        except (KeyError, IndexError):
            warn("time coordinate not found, use `coords=arange`", level=1)
            index = None

        ctrapz = cumtrapz(series, index=index, **kwargs)
        return xr.DataArray(
            ctrapz,
            dims=series.dims,
            coords=series.coords,
            name=keys._join(series.name, keys.cumtrapz),
        )

    raise TypeError("`series` not of type ndarray, Series, or DataArray?")
