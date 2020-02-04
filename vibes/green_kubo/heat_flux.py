"""compute and analyze heat fluxes"""
import numpy as np
import xarray as xr
from ase import units
from scipy.integrate import cumtrapz
from scipy.optimize import curve_fit

from vibes import dimensions as dims
from vibes import keys
from vibes.correlation import get_correlation

from .utils import Timer, talk


def gk_prefactor(volume, temperature, verbose=True):
    """convert eV/AA^2/fs to W/mK

    Args:
        volume (float): volume of the supercell in AA^3
        temperature (float): avg. temp. in K (trajectory.temperatures.mean())

    Returns:
        V / (k_B * T^2) * 1602
    """
    V = float(volume)
    T = float(temperature)
    prefactor = 1 / units.kB / T ** 2 * 1602 * V / 1000  # / 3  # * 1000
    msg = [
        f"Compute Prefactor:",
        f".. Volume:       {V:10.2f}  AA^3",
        f".. Temperature:  {T:10.2f}  K",
        f"-> Prefactor:    {prefactor:10.2f}  W/mK / (eV/AA^/ps)",
    ]
    talk(msg, verbose=verbose)
    return float(prefactor)


def _get_kappa_array(jcorr, delta="auto", verbose=True):
    """integrate HFACF to obtain kappa as xr.Dataarray and determinde avalanche time

    Args:
        jcorr: HFACF as pd.Series or xr.DataArray
        delta: mode to compute avalanche time
    Returns:
        cumulative kappa as xr.DataArray
    """
    # integrate
    talk(f"Integrate heat flux autocorrelation function cumulatively", verbose=verbose)
    talk(f".. Integrator:   `scipy.integrate.cumtrapz`", verbose=verbose)

    kappa = cumtrapz(jcorr, jcorr.time, axis=0, initial=0)

    # estimate avalanche time
    jc = jcorr[:, 0, 0] + jcorr[:, 1, 1] + jcorr[:, 2, 2]

    tmax = t_avalanche(jc.to_series(), delta=delta)

    attrs = jcorr.attrs
    attrs.update({"t_avalanche": tmax})

    # fmt: off
    da = xr.DataArray(
        kappa,
        dims=jcorr.dims,
        coords=jcorr.coords,
        attrs=attrs,
        name="kappa",
    )
    # fmt: on

    return da


def get_kappa(dataset, full=False, delta="auto", verbose=True):
    """compute heat flux autocorrelation and cumulative kappa from heat_flux_dataset

    Args:
        dataset (xr.Dataset): contains heat flux per atom
        full (bool): return correlation function per atom
        delta: compute mode for avalanche time
    Returns:
        dataset containing
            Jcorr: heat flux autocorrelation function (HFACF)
            Jcorr_II: sum_I term of HFACF
            kappa: cumulative kappa from integrating HFACF
            kappa_II: cumulative kappa from integrating Jcorr_II
    """

    if full:
        J_corr = get_heat_flux_aurocorrelation(dataset, verbose=verbose)
    else:
        flux = dataset.heat_flux.sum(axis=1)
        Nt = len(flux)
        temp = dataset.temperature.mean()
        vol = dataset.attrs["volume"]

        prefactor = gk_prefactor(vol, temp)

        flux_corr = np.zeros([Nt, 3, 3])
        for aa in range(3):
            for bb in range(3):
                corr = get_correlation(flux[:, aa], flux[:, bb])
                flux_corr[:, aa, bb] = corr

        J_corr = xr.DataArray(
            flux_corr * prefactor,
            dims=dims.stress,
            coords=flux.coords,
            attrs={"gk_prefactor": prefactor},
            name="Jcorr",
        )

    J_corr_II = get_heat_flux_aurocorrelation(dataset, verbose=verbose)
    J_corr_II.name = J_corr.name + "_II"

    kappa = _get_kappa_array(J_corr, delta=delta, verbose=verbose)
    kappa_II = _get_kappa_array(J_corr_II.sum(axis=1), delta=delta, verbose=verbose)
    kappa_II.name = kappa.name + "_II"

    datasets = (J_corr, J_corr_II, kappa, kappa_II)
    dct = {da.name: da for da in datasets}
    coords = J_corr.coords
    attrs = J_corr.attrs

    return xr.Dataset(dct, coords=coords, attrs=attrs)


def get_heat_flux_aurocorrelation(dataset, full_tensor=False, verbose=True):
    """compute heat flux autocorrelation function from heat flux Dataset

    Args:
        dataset (xarray.Dataset): the heat flux Dataset
        full_tensor (bool): return [N_t, N_a, N_a, 3, 3], otherwise [N_t, N_a, 3, 3]
    Returns:
        heat_flux_autocorrelation (xarray.DataArray [N_t, N_a, (N_a, 3,) 3]) in W/mK/fs
    """

    flux = dataset.heat_flux

    timer = Timer("Get heat flux autocorrelation from heat flux", verbose=verbose)

    if full_tensor:
        talk(".. return full tensor [N_t, N_a, N_a, 3, 3]", verbose=verbose)

    Nt, Na, _ = flux.shape
    temp = dataset.temperature.mean()
    vol = dataset.attrs["volume"]

    prefactor = gk_prefactor(vol, temp)

    if full_tensor:
        flux_corr = np.zeros([Nt, Na, Na, 3, 3])
        for II in flux.I:
            print(int(II))
            for JJ in flux.I:
                J_I = flux[:, II]
                J_J = flux[:, JJ]

                for aa in range(3):
                    for bb in range(3):
                        corr = get_correlation(J_I[:, aa], J_J[:, bb])
                        flux_corr[:, II, JJ, aa, bb] = corr
        _dims = dims.kappa
    else:
        flux_corr = np.zeros([Nt, Na, 3, 3])
        for atom in flux.I:
            J_atom = flux[:, atom]

            for aa in range(3):
                for bb in range(3):
                    corr = get_correlation(J_atom[:, aa], J_atom[:, bb])
                    flux_corr[:, atom, aa, bb] = corr
        _dims = flux.dims + (dims.kappa[-1],)

    df_corr = xr.DataArray(
        flux_corr * prefactor,
        dims=_dims,
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


def F_avalanche(series, delta="auto", verbose=True):
    """Compute Avalanche Function (windowed noise/signal ratio)

    as defined in J. Chen et al. / Phys. Lett. A 374 (2010) 2392
    See also: Parzen, Modern Probability Theory and it's Applications, Chp. 8.6, p. 378f

    Args:
        series (pandas.Series): some time resolved data series
        delta (int): no. of time steps for windowing, or `auto`
        verbose (bool): be verbose

    Returns:
        F(t, delta) = abs( sigma(series) / E(series)),
        where sigma is the standard deviation of the time series in an interval
        delta around t, and E is the expectation value around t.

    When `delta='auto'`, estimate the correlation time and use that size for binning
    the time steps for estimating std/E
    """
    if delta == "auto":
        # estimate correlation time

        exp = lambda x, y0, tau: y0 * np.exp(-x / tau)

        # where is jc drops below 0 the first time
        idx = series[series < series.max() / np.e].index[0]
        x = series[:idx]
        y0 = x.iloc[0]
        tau = idx
        bounds = ([0.75 * y0, 1.5 * y0], [0.5 * idx, 2 * idx])
        (y0, tau), _ = curve_fit(exp, x.index, x, bounds=bounds)

        delta = len(series[series.index < tau])

        if verbose:
            talk(f"Pre-est. correlation time (drop below 1/e): {idx:10.2f} fs")
            talk(f".. use {len(x)} data points to fit exponential")
            talk(f".. estimated correlation time:              {tau:10.2f} fs")
            talk(f"-> choose delta of size: {delta:20d} data points")

        assert 0.5 < idx / tau < 2, (idx, tau)

    sigma = series.rolling(window=delta, min_periods=0).std()
    E = series.rolling(window=delta, min_periods=0).mean()

    F = (sigma / E).abs().dropna()

    return F


def t_avalanche(series, Fmax=1, verbose=True, **kwargs):
    """get avalanche time for series from F_avalanche

    Args:
        Fmax (float): max. allowed value for avalanche function
    """
    F = F_avalanche(series, verbose=verbose, **kwargs)

    tmax = F[F > Fmax].index[0]

    if verbose:
        talk(f"-> avalanche time with max. F of {Fmax:2d}: {tmax:17.2f} fs")

    return tmax
