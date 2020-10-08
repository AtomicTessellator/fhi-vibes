"""compute and analyze heat fluxes"""
import scipy.signal as sl

from vibes.correlation import get_autocorrelationNd
from vibes.fourier import get_fourier_transformed
from vibes.integrate import trapz

from . import Timer, _talk


def get_velocity_autocorrelation(velocities=None, trajectory=None, verbose=True):
    """LEGACY: compute velocity autocorrelation function from xarray"""
    return get_autocorrelationNd(velocities, normalize=True, hann=False)


def get_vdos(velocities=None, hann=False, normalize=False, npad=10000, verbose=True):
    r"""compute vibrational DOS for trajectory

    vdos(w) = FT{\sum_i corr(v_i, v_i)(t)}(w)

    Args:
        velocities (xarray.DataArray [N_t, N_a, 3]): the velocities
        hann: use Hann window when computing the autocorrelation
        normalize: normalize VDOS to 1
        npad: number of zeros for zero padding
    Returns:
        vdos (xarray.DataArray [N_t, N_a, 3])
    """
    timer = Timer("Get VDOS", verbose=verbose)
    v_corr = get_autocorrelationNd(velocities, normalize=True, hann=hann)
    df_vdos = get_fourier_transformed(v_corr, npad=npad)

    if normalize:
        norm = trapz(df_vdos)
        _talk(f"Normalize with {norm}")
        df_vdos /= norm

    timer()

    return df_vdos


def simple_plot(series, file="vdos.pdf", height=-1, max_frequency=30.0, logy=False):
    """simple plot of VDOS for overview purpose

    Args:
        series (pandas.Series): Intensity vs. omega
        file (str): file to store the plot to
        height (float): minimal height to detect peaks
        max_frequency (float): max. frequency in THz
        logy (bool): use semilogy
    """
    # normalize peaks
    series /= series.max()

    # find peaks:
    if height and height > 0:
        peaks, *_ = sl.find_peaks(series, height=height)
        high_freq = series.index[peaks[-1]]
        _talk(f".. highest peak found at   {high_freq:.2f} THz")
    else:
        high_freq = max_frequency
    ax = series.plot()
    if logy:
        ax.set_yscale("log")
    ax.set_xlim([0, 1.2 * high_freq])
    ax.set_xlabel("Omega [THz]")
    fig = ax.get_figure()
    fig.savefig(file, bbox_inches="tight")
    _talk(f".. max. frequency for plot:  {high_freq:.2f} THz")
    _talk(f".. plot saved to {file}")
