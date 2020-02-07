import numpy as np
import pandas as pd
import pytest
import xarray as xr

from vibes import dimensions as dims
from vibes.correlation import get_autocorrelation, get_autocorrelationNd

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
    # test [Nt, 3] array
    a = np.zeros([5, 3])
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

    # test [Nt, Na, 3] array
    Nt, Na, Nx = 5, 4, 3
    a = np.zeros([Nt, Na, Nx])
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


if __name__ == "__main__":
    for array in arrays:
        test_type(array)
        test_autocorrelationNd()
