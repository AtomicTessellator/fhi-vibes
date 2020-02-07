import numpy as np
import pandas as pd
import pytest
import xarray as xr

from vibes import dimensions as dims
from vibes.correlation import get_autocorrelation, get_autocorrelationNd
from vibes.helpers import xtrace

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
    assert x.shape == c.shape
    assert c.dims == dims.time_vec

    c = get_autocorrelationNd(x, off_diagonal_coords=True)
    assert c.shape == (*x.shape, x.shape[-1])
    assert c.dims == dims.time_tensor

    try:
        c = get_autocorrelationNd(x, off_diagonal_atoms=True)
    except ValueError:
        pass

    # test if `off_diagonal_coords` works
    c1 = get_autocorrelationNd(x).sum(axis=-1)
    c2 = xtrace(get_autocorrelationNd(x, off_diagonal_coords=True))

    assert np.allclose(c1, c2)

    # test [Nt, Na, 3] array
    a = np.random.rand(Nt, Na, Nx)
    x = xr.DataArray(a, dims=dims.time_atom_vec)

    c = get_autocorrelationNd(x)
    assert c.shape == x.shape
    assert c.dims == dims.time_atom_vec

    c = get_autocorrelationNd(x, off_diagonal_coords=True)
    assert c.shape == (*x.shape, x.shape[-1])
    assert c.dims == dims.time_atom_tensor

    c = get_autocorrelationNd(x, off_diagonal_atoms=True)
    assert c.shape == (Nt, Na, Na, Nx, Nx)
    assert c.dims == dims.time_atom_atom_tensor

    # test if `off_diagonal_coords` works
    c1 = get_autocorrelationNd(x).sum(axis=-1)
    c2 = xtrace(get_autocorrelationNd(x, off_diagonal_coords=True))
    c3 = xtrace(get_autocorrelationNd(x, off_diagonal_atoms=True), axis1=1, axis2=2)

    assert np.allclose(c1, c2)
    assert np.allclose(c1.sum(axis=-1), xtrace(c3))


if __name__ == "__main__":
    for array in arrays:
        test_type(array)
        test_autocorrelationNd()
