"""test fourier transform"""

import numpy as np
import pandas as pd
from scipy import signal as sl

from vibes.fourier import get_fourier_transformed

# create signal with 100 THz
omega = 100
to_THZ = np.pi * 2 / 1000

t = np.linspace(0, 1000, 10000)
y = np.sin(omega * to_THZ * t)  # + 10

ndarray = y
series = pd.Series(data=y, index=t)


def test_fourier(f=series, t=t, omega=omega):
    """Check if the fourier transform has a peak where expected"""
    # get_fourier only accepts xarray.DataArrays
    ft = get_fourier_transformed(f.to_xarray()).real.to_series()

    peak_index = sl.find_peaks(ft)[0]
    peak_omega = ft.index[peak_index]

    ft.to_csv("test.csv")

    assert np.allclose(peak_omega, omega, rtol=1e-3)


if __name__ == "__main__":
    test_fourier()
