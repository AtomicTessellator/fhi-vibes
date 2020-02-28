"""test fourier transform"""
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from scipy import signal as sl

from vibes.fourier import get_fourier_transformed

# create signal with 100 THz
omega = 100
to_THZ = np.pi * 2 / 1000

t = np.linspace(0, 1000, 10000)
y = np.sin(omega * to_THZ * t) + 10

ndarray = y
series = pd.Series(y, t)
xarray = xr.DataArray(ndarray)

arrays = (ndarray, series, xarray)


def test_fourier(f=series, t=t, omega=omega):
    """check if the fourier transform has a peak where expected"""
    ft = get_fourier_transformed(f)

    peak_index = sl.find_peaks(ft)[0]
    peak_omega = ft.index[peak_index]

    print(peak_omega / omega)

    assert np.allclose(peak_omega, omega, rtol=1e-3)


@pytest.mark.parametrize("array", arrays)
def test_type(array):
    fft = get_fourier_transformed(array)
    assert type(fft) == type(array)
    # assert len(fft) == len(array), (len(fft), len(array))


if __name__ == "__main__":
    for array in arrays:
        test_type(array)
    test_fourier()
