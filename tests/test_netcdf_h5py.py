# https://github.com/Unidata/netcdf4-python/issues/694
from pathlib import Path

import netCDF4

try:
    import h5py  # noqa: F401
except ModuleNotFoundError:
    pass

parent = Path(__file__).parent


def test_h5py():
    try:
        netCDF4.Dataset(parent / "trajectory" / "trajectory.nc", "r")
    except OSError:
        msg = (
            "Please check h5py installation, see: \n",
            "https://github.com/Unidata/netcdf4-python/issues/694",
        )
        raise RuntimeError(msg) from None
