"""test fourier transform"""
import numpy as np
import pandas as pd
from scipy import signal as sl

from vibes.fourier import get_fourier_transformed

# create signal with 100 THz
omega = 100
to_THZ = np.pi * 2 / 1000

t = np.linspace(0, 1000, 10000)
f = pd.Series(np.sin(omega * to_THZ * t) + 10, t)


def test_fourier(f=f, t=t, omega=omega):
    """check if the fourier transform has a peak where expected"""
    ft = get_fourier_transformed(f)

    peak_index = sl.find_peaks(ft)[0]
    peak_omega = ft.omega[peak_index]

    print(peak_omega / omega)

    assert np.allclose(peak_omega, omega, rtol=1e-3)


if __name__ == "__main__":
    test_fourier()
