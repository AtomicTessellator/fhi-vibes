"""Fourier Transforms"""
import numpy as np
import pandas as pd
import xarray as xr

from vibes import dimensions as dims
from vibes import keys
from vibes.helpers import Timer, talk, warn
from vibes.konstanten.einheiten import THz_to_cm

_prefix = "Fourier"


def get_timestep(times, tol=1e-9):
    """get time step from a time series and check for glitches"""
    d_times = np.asarray((times - np.roll(times, 1))[1:])
    timestep = np.mean(d_times)

    assert np.all(np.abs(d_times - timestep) < tol)

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
        talk(msg, prefix=_prefix)

    if to_cm:
        w *= THz_to_cm

    return w


def get_fft(series):
    """run `np.fft.fft` on time series

    https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.hfft.html

    Args:
        series (np.ndarray [N_t, ...]): time series, first dimension is the time axis

    Returns:
        np.ndarray ([N_t, ...]): Fourier transform of series ([: N_t // 2])
    """

    velocities = np.asarray(series).copy()

    N = series.shape[0]

    velocities = velocities.swapaxes(-1, 0)
    velocities = np.fft.fft(velocities, axis=-1)

    return velocities.swapaxes(-1, 0)[: N // 2]


def get_fourier_transformed(series, verbose=True):
    """Perform Fourier Transformation of Series/DataArray

    Args:
        series ([N_t, ...]): pandas.Series/xarray.DataArray with `time` axis in fs
        verbose (bool): be verbose
    Return:
        DataArray ([N_t, ...]): FT(series) with `omega` axis in THz
    """
    timer = Timer("Compute FFT", verbose=verbose, prefix=_prefix)

    fft = get_fft(series)

    if isinstance(series, np.ndarray):
        result = fft
    elif isinstance(series, pd.Series):
        omegas = get_frequencies(times=series.index, verbose=verbose)
        result = pd.Series(fft, index=omegas)
        result.index.name = keys.omega
    elif isinstance(series, xr.DataArray):
        try:
            times = np.asarray(series[keys.time])
            omegas = get_frequencies(times=times, verbose=verbose)
        except (KeyError, IndexError):
            warn(f"time coordinate not found, use `coords=arange`", level=1)
            omegas = get_frequencies(times=np.arange(len(series)), verbose=verbose)

        da = xr.DataArray(
            fft,
            dims=(keys.omega, *series.dims[1:]),
            coords={keys.omega: omegas},
            name=keys._join(series.name, keys.fourier_transform),
        )
        result = da
    else:
        raise TypeError("`series` not of type ndarray, Series, or DataArray?")

    timer()

    return result


def get_power_spectrum(
    flux: xr.DataArray = None, prefactor: float = 1.0, verbose: bool = True
) -> xr.DataArray:
    """compute power spectrum for given flux

    Args:
        flux (xr.DataArray, optional): flux [Nt, 3]. Defaults to xr.DataArray.
        prefactor (float, optional): prefactor
        verbose (bool, optional): be verbose

    Returns:
        xr.DataArray: heat_flux_power_spectrum
    """
    kw = {"verbose": verbose}
    timer = Timer("Get power spectrum from flux", **kw)

    Jw = get_fourier_transformed(flux.dropna(dims.time), **kw)

    Jspec = abs(Jw)

    timer()

    return Jspec
