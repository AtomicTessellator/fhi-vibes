"""compute and analyze heat fluxes"""
import numpy as np
import scipy.signal as sl
from scipy.integrate import cumtrapz
import xarray as xr

from ase import units
# from hilde.fourier import compute_sed, get_frequencies, get_timestep
from . import Timer, talk


def gk_prefactor(volume, temperature):
    """convert eV/AA^2/ps to W/mK

    Args:
        volume (float): volume of the supercell in AA^3
        temperature (float): avg. temp. in K (trajectory.temperatures.mean())

    Returns:
        V / (3 * k_B * T^2) * 1602
    """
    V = float(volume)
    T = float(temperature)
    prefactor = 1 / units.kB / T ** 2 * 1602 * V / 3  # * 1000
    talk(
        [
            f"Compute Prefactor:",
            f".. Volume:       {V:10.2f}  AA^3",
            f".. Temperature:  {T:10.2f}  K",
            f"-> Prefactor:    {prefactor:10.2f}  W/mK / (eV/AA^/ps)",
        ]
    )
    return float(prefactor)


def get_heat_flux_aurocorrelation(dataset, verbose=True):
    """compute heat flux autocorrelation function from heat flux Dataset

    Args:
        dataset (xarray.Dataset): the heat flux Dataset
    Returns:
        heat_flux_autocorrelation (xarray.DataArray [N_t, N_a, 3]) in W/mK/fs
    """

    flux = dataset.heat_flux

    timer = Timer("Get heat flux autocorrelation from heat flux", verbose=verbose)

    Nt = len(flux.time)
    temp = dataset.temperature.mean()
    vol = dataset.attrs["volume"]

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
        attrs={"gk_prefactor": prefactor},
        name="heat flux autocorrelation",
    )

    timer()

    return df_corr


def get_cumulative_kappa(dataset, analytic=False, verbose=True):
    """compute cumulative kappa(T) with T the integration boundary in fs

    Args:
        dataset (xarray.Dataset): the heat flux Dataset with time in fs
        analytic (bool): integrate analytially by interpolating the flux

    Returns:
        cumulative_kappa (xarray.DataArray [N_t, N_a, 3]) in W/mK
    """

    # J_corr in .../ps
    J_corr = get_heat_flux_aurocorrelation(dataset)

    # dt in ps
    times = dataset.coords["time"] / 1000
    # dt = get_timestep(times)
    # d_times = (times - np.roll(times, 1))[1:]

    # integrate
    talk(f"Integrate heat flux autocorrelation function cumulatively")
    talk(f".. Integrator:   `scipy.integrate.cumtrapz`")
    talk(f".. analytic:      {analytic}")

    kappa = cumtrapz(J_corr, times, axis=0, initial=0)

    # fmt: off
    da = xr.DataArray(
        kappa,
        dims=J_corr.dims,
        coords=J_corr.coords,
        attrs=J_corr.attrs,
        name="kappa",
    )
    # fmt: on

    return da


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
