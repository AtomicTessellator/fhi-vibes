import numpy as np
import xarray as xr


def xtrace(array: xr.DataArray, axis1: int = -2, axis2: int = -1) -> xr.DataArray:
    """apply `np.trace` to xarray.DataArray

    Args:
        array (xr.DataArray): [description]
        axis1 (int, optional): [description]. Defaults to -2.
        axis2 (int, optional): [description]. Defaults to -1.

    Returns:
        xr.DataArray: DataArray with trace applied
    """
    data = np.trace(array, axis1=axis1, axis2=axis2)

    da = xr.DataArray(data, dims=array.dims[:-2], coords=array.coords)

    return da
