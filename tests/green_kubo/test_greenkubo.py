""" test green kubo cumulative kappa"""
from pathlib import Path

import xarray as xr

from vibes import green_kubo as gk
from vibes import keys


parent = Path(__file__).parent

ds = xr.load_dataset(parent / "test.nc")


def test_get_hf_data():
    hfacf, kappa = gk.get_hf_data(ds[keys.heat_flux])

    for array in (hfacf, kappa):
        assert isinstance(array, xr.DataArray)


def test_get_filtered():
    array = ds[keys.velocities]

    array_filtered = gk.get_filtered(array, window=2)

    assert isinstance(array_filtered, xr.DataArray)
    assert array.shape == array_filtered.shape

    # check antisymmetric
    array_filtered = gk.get_filtered(array, window=2, antisymmetric=True)

    assert isinstance(array_filtered, xr.DataArray)
    assert array.shape == array_filtered.shape


def test_get_gk_dataset():
    gk.get_gk_dataset(ds, filter_threshold=0.8)
