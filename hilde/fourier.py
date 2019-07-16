"""Fourier Transforms"""
import numpy as np
from hilde.konstanten.einheiten import THz_to_cm
from hilde.helpers import talk


def get_timestep(times):
    """get time step from a time series"""
    d_times = (times - np.roll(times, 1))[1:]
    timestep = np.mean(d_times)

    if any(d_times - timestep > 1e-9):
        msg = f"delta times - timestep: {d_times - timestep}\n"
        msg += f"timesteps uneven? Inspect times!"
        raise ValueError(msg)

    return timestep


def get_frequencies(
    N=None, dt=None, times=None, fs_factor=1, to_cm=False, verbose=False
):
    """compute the frequency domain in THz for signal with fs resolution

    Args:
        N (int): Number of time steps
        dt (float): time step in fs
        times (ndarray): time series in fs (or converted to fs via `fs_factor`)
        fs_factor (float): convert timestep to fs by `dt / fs_factor` (for ps: 1000)
        to_cm (bool): return in inverse cm instead
        verbose (bool): be informative

    Returns:
        frequencies in THz (ndarray)
    """
    if N and dt:
        dt = dt / fs_factor
    elif times is not None:
        N = len(times)
        dt = get_timestep(times) / fs_factor

    # Frequencies in PetaHz
    max_freq = 1.0 / (2.0 * dt)
    w = np.linspace(0.0, max_freq, N // 2)
    # Frequencies in THz
    w *= 1000
    dw = w[1] - w[0]

    if verbose:
        msg = f".. get frequencies for time series\n"
        msg += f".. timestep:               {np.asarray(dt)} fs\n"
        msg += f"-> maximum frequency:      {np.max(w):.5} THz\n"
        msg += f".. Number of steps:        {N}\n"
        msg += f"-> frequency resolution:   {dw:.5f} THz\n"
        talk(msg)

    if to_cm:
        w *= THz_to_cm

    return w


def compute_sed(series):
    """ Computes spectral density for a time series

    Args:
        series (np.ndarray [N_t, ...]): time series, first dimension is the time axis

    Returns:
        Fourier transform of series
    """

    velocities = series.copy()

    N = series.shape[0]

    velocities = velocities.swapaxes(-1, 0).copy()
    velocities = np.fft.fft(velocities, axis=-1)

    return velocities.swapaxes(-1, 0)[: N // 2]
