"""compute and analyze heat fluxes"""
import numpy as np
import scipy.signal as sl
import xarray as xr

from ase import units
from hilde.fourier import compute_sed, get_frequencies, get_timestep
from hilde.helpers import progressbar, talk, Timer


def gk_prefactor(volume, temperature):
    """convert eV/AA^2/ps to W/mK

    Args:
        volume (float): volume of the supercell in AA^3
        temperature (float): avg. temp. in K (trajectory.temperatures.mean())

    Returns:
        V / (3 * k_B * T^2) * 1602
    """
    V = volume
    T = temperature
    prefactor = 1 / units.kB / T ** 2 * 1602 * V / 3
    return prefactor


def atoms_with_stresses(trajectory):
    """return new trajectory with atoms that have stresses computed"""
    return trajectory.with_stresses


def average_atomic_stress(trajectory, verbose=True):
    """compute average atomic stress from trajectory"""
    timer = Timer("Compute average atomic stress:", verbose=verbose)

    atomic_stresses = []
    for a in progressbar(trajectory):
        V = a.get_volume()

        atomic_stress = a.get_stresses() / V
        atomic_stresses.append(atomic_stress)

    timer()

    return np.mean(atomic_stresses, axis=0)


def get_heat_flux(trajectory, return_avg=False):
    """compute heat fluxes from TRAJECTORY and return as xarray

    Args:
        trajectory: list of atoms objects
        return_avg (bool): return flux per atom

    Returns:
        flux ([N_t, N_a, 3]): the time resolved heat flux in eV/AA**3/ps
        avg_flux (optional, [N_t, N_a, 3]): high frequency part of heat flux
    """

    # get atoms that have stresses computed
    atoms_w_stress = atoms_with_stresses(trajectory)

    # 1) compute average stresses
    avg_stresses = average_atomic_stress(atoms_w_stress)

    # 2) compute J_avg from average stresses and dJ from variance
    timer = Timer("Compute heat fluxes:")
    fluxes = []
    avg_fluxes = []
    times = []
    for a in progressbar(atoms_w_stress):
        times.append(a.info["dt_fs"] * a.info["nsteps"])
        V = a.get_volume()
        stresses = a.get_stresses() / V

        ds = stresses - avg_stresses

        # velocity in \AA / ps
        vs = a.get_velocities() * units.fs * 1000

        fluxes.append((ds @ vs[:, :, None]))
        avg_fluxes.append((avg_stresses @ vs[:, :, None]))

    timestep = get_timestep(times)
    talk(f".. time step is {timestep:.2f} fs")

    times = np.array(times) - min(times)

    if return_avg:
        fluxes = np.array(avg_fluxes).squeeze()
        name = "avg_heat_flux"
    else:
        fluxes = np.array(fluxes).squeeze()
        name = "heat_flux"

    df = xr.DataArray(
        fluxes, dims=["time", "atom", "coord"], coords={"time": times}, name=name
    )

    timer()

    return df


def get_heat_flux_aurocorrelation(trajectory, verbose=True):
    """compute heat flux autocorrelation function from trajectory

    Args:
        trajectory: list of atoms objects
    Returns:
        heat_flux_autocorrelation (xarray.DataArray [N_t, N_a, 3]) in W/mK/ps
    """
    traj = trajectory.with_stresses

    flux = get_heat_flux(traj)

    timer = Timer("Get heat flux autocorrelation from trajectory", verbose=verbose)

    Nt = len(flux.time)
    temp = traj.temperatures.mean()
    vol = traj[0].get_volume()

    prefactor = gk_prefactor(vol, temp)

    flux_corr = np.zeros_like(flux)
    for atom in flux.atom:
        J_atom = flux[:, atom]

        for xx in range(3):
            corr = sl.correlate(J_atom[:, xx], J_atom[:, xx])[Nt - 1 :] / Nt
            flux_corr[:, atom, xx] = corr

    df_corr = xr.DataArray(
        flux_corr * prefactor,
        dims=flux.dims,
        coords=flux.coords,
        name="heat flux autocorrelation",
    )

    timer()

    return df_corr


def get_cumulative_kappa(trajectory, verbose=True):
    """get kappa(T)"""

    traj = trajectory.with_stresses

    # dt in ps
    dt = get_timestep(traj.times) / 1000

    # J_corr in .../ps
    J_corr = get_heat_flux_aurocorrelation(traj)

    kappa = (J_corr * dt).cumsum(axis=0)

    kappa.name = "kappa"

    return kappa


def F_avalanche(series, delta=20):
    """Compute Avalanche Function (windowed noise/signal ratio)

    as defined in J. Chen et al. / Phys. Lett. A 374 (2010) 2392
    See also: Parzen, Modern Probability Theory and it's Applications, Chp. 8.6, p. 378f

    Args:
        series (pandas.Series): some time resolved data series
        delta (int): no. of time steps for windowing

    Returns:
        F(t, delta) = abs( sigma(series) / E(series)),
        where sigma is the standard deviation of the time series in an interval
        delta around t, and E is the expectation value around t.
    """
    sigma = series.rolling(window=delta, min_periods=0).std()
    E = series.rolling(window=delta, min_periods=0).mean()

    F = (sigma / E).abs()

    return F
