"""wrappers for scipy.integrate functions"""

import numpy as np
import xarray as xr
from scipy import integrate as si

from vibes import keys
from vibes.helpers import Timer  # , talk

_prefix = "Integration"
Timer.prefix = _prefix


def _cumtrapz(series, index=None, axis=0, initial=0):
    """wrap `scipy.integrate.cumtrapz`"""
    array = np.asarray(series)
    ct = si.cumtrapz(array, x=index, axis=axis, initial=initial)
    return ct


def get_cumtrapz(series, **kwargs):
    """Compute cumulative trapezoid integral of Series/DataArray

    Return:
        DataArray
    """
    xarray = True
    try:
        index = series[keys.time]
    except KeyError:
        xarray = False
        index = series.index

    ct = _cumtrapz(series, index=index, **kwargs)

    # fmt: off
    da = xr.DataArray(
        ct,
        dims=[keys.time],
        coords={keys.time: index}, name=keys.cumtrapz,
    )
    # fmt: on

    if not xarray:
        return da.to_series()
    return da
