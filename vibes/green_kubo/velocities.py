"""compute and analyze heat fluxes"""
import scipy.signal as sl
import xarray as xr

from vibes.correlation import get_autocorrelationNd
from vibes.fourier import get_fft, get_frequencies
from vibes.helpers import Timer, talk
from vibes.trajectory.dataset import get_velocities_dataarray


def get_velocity_autocorrelation(velocities=None, trajectory=None, verbose=True):
    """LEGACY: compute velocity autocorrelation function from xarray"""
    return get_autocorrelationNd(velocities, normalize=True, window=False)


def get_vdos(velocities=None, trajectory=None, verbose=True):
    r"""compute vibrational DOS for trajectory

    vdos(w) = FT{\sum_i corr(v_i, v_i)(t)}(w)

    Args:
        velocities (xarray.DataArray [N_t, N_a, 3]): the velocities
        trajectory: list of atoms objects
    Returns:
        vdos (xarray.DataArray [N_t, N_a, 3])
    """
    if velocities is None and trajectory is not None:
        velocities = get_velocities_dataarray(trajectory, verbose=verbose)

    v_corr = get_autocorrelationNd(velocities, normalize=True, window=False)

    timer = Timer("Get VDOS", verbose=verbose)

    omegas = get_frequencies(times=v_corr.time, verbose=verbose)

    v_spec = get_fft(v_corr.data)

    # fmt: off
    df_vdos = xr.DataArray(
        v_spec,
        dims=["omega", *v_corr.dims[1:]],
        coords={"omega": omegas}, name="vdos",
    )
    # fmt: on

    timer()

    return df_vdos


def simple_plot(series, file="vdos.pdf", height=-1, max_frequency=30.0):
    """simple plot of VDOS for overview purpose

    Args:
        series (pandas.Series): Intensity vs. omega
        file (str): file to store the plot to
        height (float): minimal height to detect peaks
        max_frequency (float): max. frequency in THz
    """
    # normalize peaks
    series /= series.max()

    # find peaks:
    if height and height > 0:
        peaks, *_ = sl.find_peaks(series, height=height)
        high_freq = series.index[peaks[-1]]
        talk(f".. highest peak found at   {high_freq:.2f} THz")
    else:
        high_freq = max_frequency
    ax = series.plot()
    ax.set_xlim([0, 1.2 * high_freq])
    ax.set_xlabel("Omega [THz]")
    fig = ax.get_figure()
    fig.savefig(file, bbox_inches="tight")
    talk(f".. max. frequency for plot:  {high_freq:.2f} THz")
    talk(f".. plot saved to {file}")
