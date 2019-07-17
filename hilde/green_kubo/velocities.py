"""compute and analyze heat fluxes"""
import numpy as np
import scipy.signal as sl
import xarray as xr

from hilde.trajectory import reader, Trajectory
from hilde.fourier import compute_sed, get_frequencies
from hilde.helpers import Timer, talk


def get_velocities(trajectory, verbose=True):
    """extract velocties from TRAJECTORY  and return as xarray.DataArray

    Args:
        trajectory (Trajectory or Path): list of atoms objects or where to find them
    Returns:
        velocities (xarray.DataArray [N_t, N_a, 3])
    """
    timer = Timer("Get velocities from trajectory", verbose=verbose)

    if not isinstance(trajectory, Trajectory):
        trajectory = reader(trajectory)

    metadata = {
        "time unit": "fs",
        "timestep": trajectory.timestep,
        "atoms": trajectory[0].get_chemical_symbols(),
    }

    times = trajectory.times
    velocities = np.array([a.get_velocities() for a in trajectory])

    df = xr.DataArray(
        velocities,
        dims=["time", "atom", "coord"],
        coords={"time": times},
        name="velocities",
        attrs=metadata,
    )

    timer()

    return df


def get_velocity_aurocorrelation(velocities=None, trajectory=None, verbose=True):
    """compute velocity autocorrelation function from xarray

    Args:
        velocities (xarray.DataArray [N_t, N_a, 3]): the velocities
        trajectory: list of atoms objects
    Returns:
        velocity_autocorrelation (xarray.DataArray [N_t, N_a, 3])
    """
    if velocities is None and trajectory is not None:
        velocities = get_velocities(trajectory, verbose=verbose)

    timer = Timer("Get velocity autocorrelation", verbose=verbose)

    Nt = len(velocities.time)

    v_atom_corr = np.zeros_like(velocities)
    for atom in velocities.atom:
        v_atom = velocities[:, atom]

        for xx in range(3):
            corr = sl.correlate(v_atom[:, xx], v_atom[:, xx])[Nt - 1 :] / Nt
            v_atom_corr[:, atom, xx] = corr

    df_corr = xr.DataArray(
        v_atom_corr,
        dims=velocities.dims,
        coords=velocities.coords,
        name="velocity_autocorrelation",
    )

    timer()

    return df_corr


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
        velocities = get_velocities(trajectory, verbose=verbose)

    v_corr = get_velocity_aurocorrelation(velocities)

    timer = Timer("Get VDOS", verbose=verbose)

    omegas = get_frequencies(times=v_corr.time, verbose=verbose)

    v_spec = compute_sed(v_corr.data)

    # fmt: off
    df_vdos = xr.DataArray(
        v_spec,
        dims=["omega", *v_corr.dims[1:]],
        coords={"omega": omegas}, name="vdos",
    )
    # fmt: on

    timer()

    return df_vdos


def simple_plot(series, file="vdos.pdf"):
    """simple plot of VDOS for overview purpose

    Args:
        series (pandas.Series): Intensity vs. omega
    """
    # find peaks:
    peaks, *_ = sl.find_peaks(series)
    high_freq = series.index[peaks[-1]]
    ax = series.plot()
    ax.set_xlim([0, 1.2 * high_freq])
    ax.set_xlabel("Omega [THz]")
    fig = ax.get_figure()
    fig.savefig(file, bbox_inches="tight")
    talk(f".. highest peak found at   {high_freq:.2f} THz")
    talk(f".. plot saved to {file}")
