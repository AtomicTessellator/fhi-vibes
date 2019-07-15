"""compute and analyze heat fluxes"""
import numpy as np
import scipy.signal as sl
import xarray as xr

from hilde.fourier import compute_sed, get_frequencies
from hilde.helpers import Timer


def get_velocites(trajectory, verbose=True):
    """extract velocties from TRAJECTORY  and return as xarray.DataSet

    Args:
        trajectory: list of atoms objects

    Returns:
        velocities (xarray.DataArray [N_t, N_a, 3])
    """
    timer = Timer("Get velocities from trajectory", verbose=verbose)

    times = trajectory.times

    velocities = np.array([a.get_velocities() for a in trajectory])

    df = xr.DataArray(
        velocities,
        dims=["time", "atom", "coord"],
        coords={"time": times},
        name="velocities",
    )

    timer()

    return df


def get_velocity_aurocorrelation(trajectory, verbose=True):
    """compute velocity autocorrelation function from trajectory

    Args:
        trajectory: list of atoms objects
    Returns:
        velocity_autocorrelation (xarray.DataArray [N_t, N_a, 3])
    """
    velocities = get_velocites(trajectory)

    timer = Timer("Get velocity autocorrelation from trajectory", verbose=verbose)

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


def get_vdos(trajectory, verbose=True):
    r"""compute vibrational DOS for trajectory

    vdos(w) = FT{\sum_i corr(v_i, v_i)(t)}(w)

    Args:
        trajectory: list of atoms objects

    Returns:
        vdos (xarray.DataArray [N_t, N_a, 3])
    """
    v_corr = get_velocity_aurocorrelation(trajectory)

    timer = Timer("Get VDOS from trajectory", verbose=verbose)

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
