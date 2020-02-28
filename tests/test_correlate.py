import numpy as np
import pandas as pd
import pytest
import xarray as xr

from vibes import dimensions as dims
from vibes.correlation import get_autocorrelation, get_autocorrelationNd
from vibes.helpers.xarray import xtrace

np.random.seed(4)

ndarray = np.array([0.0, 1.0, 1.0]).repeat(4)
series = pd.Series(ndarray)
xarray = xr.DataArray(ndarray)

arrays = (ndarray, series, xarray)


@pytest.mark.parametrize("array", arrays)
def test_type(array):
    corr = get_autocorrelation(array)
    assert type(corr) == type(array)
    assert len(corr) == len(array)


def test_autocorrelationNd():
    Nt, Na, Nx = 5, 4, 3

    # test [Nt, 3] array
    a = np.random.rand(Nt, Nx)
    x = xr.DataArray(a, dims=dims.time_vec)

    c = get_autocorrelationNd(x)
    assert c.shape == (Nt, Nx)
    assert c.dims == dims.time_vec

    c1 = get_autocorrelationNd(x, off_diagonal=True, distribute=False)
    assert c1.shape == (Nt, Nx, Nx)

    c2 = get_autocorrelationNd(x, off_diagonal=True, distribute=False)
    assert np.allclose(c1, c2)

    # test if `off_diagonal_coords` works
    c1 = get_autocorrelationNd(x).sum(axis=-1)
    c2 = xtrace(get_autocorrelationNd(x, off_diagonal=True))

    assert np.allclose(c1, c2)

    # test [Nt, Na, 3] array
    a = np.random.rand(Nt, Na, Nx)
    x = xr.DataArray(a, dims=dims.time_atom_vec)

    c = get_autocorrelationNd(x)
    assert c.shape == x.shape
    assert c.dims == dims.time_atom_vec


if __name__ == "__main__":
    for array in arrays:
        test_type(array)
        test_autocorrelationNd()
