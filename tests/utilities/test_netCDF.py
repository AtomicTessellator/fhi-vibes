"""test netCDF capabilites"""

from pathlib import Path
import numpy as np
import xarray as xr

file = Path("test.nc")

times = np.linspace(0, 10, 10)
sine = np.sin(2 * np.pi * times)
metadata = {
    "timestep": times[1] - times[0],
    "creation_data": "19-07-17",
    "str_list": ["a", "b"],
    "int_list": [2, 3],
}

df = xr.DataArray(
    sine, dims="time", coords={"time": times}, name="sine", attrs=metadata
)


def test_write():
    """test writing an xarray.DataArray file"""
    df.to_netcdf(file)


def test_read():
    """test reading an xarray.DataArray file"""
    new_df = xr.open_dataarray(file)

    for key, val in df.attrs.items():
        assert np.all(val == new_df.attrs[key]), (df.attrs, new_df.attrs)
        assert all(df == new_df), (df, new_df)

    file.unlink()


if __name__ == "__main__":
    test_write()
    test_read()
