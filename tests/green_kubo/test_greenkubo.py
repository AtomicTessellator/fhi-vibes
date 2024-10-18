"""test green kubo cumulative kappa"""

from pathlib import Path

import xarray as xr
import numpy as np

from vibes import green_kubo as gk
from vibes import keys

parent = Path(__file__).parent

ds = xr.load_dataset(parent / "test.nc")
gk_ds = xr.load_dataset(parent / "greenkubo_test.nc")


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


def test_get_gk_dataset(ref_gk_ds=gk_ds):
    _keys = (
        keys.hf_acf,
        keys.kappa_cumulative,
    )

    GK_DS = gk.get_gk_dataset(ds)

    for key in _keys:
        assert np.allclose(ref_gk_ds.get(key).data, GK_DS.get(key).data)


if __name__ == "__main__":
    test_get_hf_data()
    test_get_filtered()
    test_get_gk_dataset()
