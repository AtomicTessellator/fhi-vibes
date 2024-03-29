import numpy as np
import pandas as pd
import pytest
import xarray as xr

from vibes.integrate import get_cumtrapz

ndarray = np.array([0.0, 1.0, 1.0]).repeat(4)
series = pd.Series(ndarray)
xarray = xr.DataArray(ndarray)

arrays = (ndarray, series, xarray)


@pytest.mark.parametrize("array", arrays)
def test_type(array):
    corr = get_cumtrapz(array)
    assert isinstance(corr, type(array))
    assert len(corr) == len(array)


if __name__ == "__main__":
    for array in arrays:
        test_type(array)
