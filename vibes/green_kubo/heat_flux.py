"""compute and analyze heat fluxes"""
import numpy as np
import xarray as xr

# import xrcache as xc
from ase import units

from vibes import defaults
from vibes import dimensions as dims
from vibes import keys
from vibes.correlation import get_autocorrelationNd
from vibes.fourier import get_fourier_transformed, get_power_spectrum
from vibes.helpers import warn
from vibes.integrate import get_cumtrapz

from .utils import Timer
from .utils import _talk as talk


def gk_prefactor(
    volume: float, temperature: float, fs_factor: float = 1, verbose: bool = True
) -> float:
    """convert eV/AA^2/fs to W/mK

    Args:
        volume (float): volume of the supercell in AA^3
        temperature (float): avg. temp. in K (trajectory.temperatures.mean())
        fs_factor (float): time * fs_factor = time in fs

    Returns:
        V / (k_B * T^2) * 1602
    """
    V = float(volume)
    T = float(temperature)
    prefactor = 1 / units.kB / T ** 2 * 1.602 * V / fs_factor  # / 1000
    msg = [
        f"Compute Prefactor:",
        f".. Volume:        {V:10.2f}  AA^3",
        f".. Temperature:   {T:10.2f}  K",
        f".. factor to fs.: {fs_factor:10.5f}",
        f"-> Prefactor:     {prefactor:10.2f}  W/mK / (eV/AA^/fs)",
    ]
    talk(msg, verbose=verbose)
    return float(prefactor)


def get_gk_prefactor_from_dataset(dataset: xr.Dataset, verbose: bool = True) -> float:
    """get the GK prefactor for the dataset, wraps `gk_prefactor`"""
    volume = dataset.attrs[keys.volume]
    temperature = dataset[keys.temperature].mean()
    return gk_prefactor(volume=volume, temperature=temperature, verbose=verbose)


# @xc.stored(verbose=True)
def get_kappa_cumulative_dataset(
    array: xr.Dataset,
    delta: str = "auto",
    verbose: bool = True,
    discard: int = 0,
    Fmax: int = defaults.F_max,
    F_window: int = defaults.F_window,
    **kw_correlate,
) -> xr.Dataset:
    """compute heat flux autocorrelation and cumulative kappa from heat_flux_dataset

    Args:
        dataset (xr.Dataset): contains heat flux per atom
        delta: compute mode for avalanche time
        discard: discard this many time steps
        Fmax: max. value for avalanche function
        F_window: window size for avalanche function
        kw_correlate (dict): kwargs for `correlate`
    Returns:
        dataset containing
            hfacf, hfacf_scalar, kappa_cumulative, kappa_cumulative_scalar
    """
    dataset = array.copy()
    pref = get_gk_prefactor_from_dataset(dataset, verbose=verbose)

    kw = {"prefactor": pref, "verbose": verbose, **kw_correlate}

    def get_jcorr(da, avg=True, discard=discard, cache=False):
        """local heper function to compute hfacf"""
        flux = da.dropna(dims.time)[discard:]
        avg_flux = flux.mean(axis=0)
        if avg:
            return get_heat_flux_aurocorrelation(flux - avg_flux, **kw)
        return get_heat_flux_aurocorrelation(flux, **kw)

    # compute hfacf
    JJ = dataset[keys.heat_flux].dropna(keys.time)
    J_corr = get_jcorr(JJ)
    # compute cumulative kappa
    kappa = get_cumtrapz(J_corr)
    dataarrays = [JJ, J_corr, kappa]

    if keys.heat_flux_aux in dataset:
        talk(".. aux. heat flux found in dataset")
        J_aux = dataset[keys.heat_flux_aux].dropna(keys.time)
        J_total = JJ + J_aux
        J_total.name = keys.heat_flux_total
        J_corr_total = get_jcorr(J_total)
        kappa_aux = get_cumtrapz(J_corr_total)
        dataarrays += [J_total, J_corr_total, kappa_aux]

    # create dataset
    dct = {da.name: da for da in dataarrays}
    coords = J_corr.coords
    attrs = J_corr.attrs

    dct.update({keys.gk_prefactor: pref})
    dct.update(dataset[[keys.volume, keys.temperature]])

    # DEV: compute separately for each component
    talk(f"Compute avalanche function with F_max={Fmax} and F_window={F_window}")

    # compute avalanche function F
    kw = {"min_periods": 5, "center": True}
    R = J_corr.rolling({keys.time: F_window}, **kw)
    E = R.mean()
    S = R.std()
    F = abs(S / E).dropna(keys.time)

    # get avalanche times per component and save the respective kappa
    ts = np.zeros([3, 3])
    ks = np.zeros([3, 3])

    for (ii, jj) in np.ndindex(ts.shape):
        f = F[:, ii, jj]
        k = kappa[:, ii, jj]
        try:
            t = f[f > Fmax].time[0]
        except IndexError:
            warn("Avalanche time not found", level=1)
            t = f.time[-1]
        ts[ii, jj] = t
        ks[ii, jj] = k.sel({keys.time: t})

    dct.update({keys.avalanche_function: F})
    attrs.update({"F_max": Fmax, "F_window": F_window})

    dct.update({keys.kappa: (dims.tensor, ks)})
    dct.update({keys.time_avalanche: (dims.tensor, ts)})

    # report
    if verbose:
        k_diag = np.diag(ks)
        talk(["Avalanche times (fs):", *np.array2string(ts, precision=3).split("\n")])
        talk(f"Kappa is:       {np.mean(k_diag):.3f} +/- {np.std(k_diag) / 3**.5:.3f}")
        talk(["Kappa^ab is: ", *np.array2string(ks, precision=3).split("\n")])

    # power spectrum
    Jw1 = get_fourier_transformed(J_corr).real
    Jw2 = get_fourier_transformed(J_corr_total, verbose=False).real
    dct.update({keys.heat_flux_power_spectrum: Jw1})
    dct.update({keys.heat_flux_total_power_spectrum: Jw2})

    DS = xr.Dataset(dct, coords=coords, attrs=attrs)

    return DS


def get_heat_flux_aurocorrelation(
    array: xr.DataArray,
    prefactor: float = 1.0,
    verbose: bool = True,
    assert_vanishing_mean: bool = False,
    **kw_correlate,
) -> xr.DataArray:
    """compute heat flux autocorrelation function from heat flux Dataset

    Args:
        flux (xr.DataArray, optional): heat flux [Nt, 3. Defaults to xr.DataArray.
        fluxes (xr.DataArray, optional): atomic heat flux [Nt, Na, 3]
        prefactor (float, optional): prefactor to convert to W/mK/fs
        verbose (bool, optional): be verbose
        assert_vanishin_mean (bool): Assert that time average is vanishing
        kw_correlate (dict): kwargs for `correlate`

    Returns:
        xr.DataArray: heat_flux_autocorrelation in W/mK/fs
    """
    timer = Timer("Get heat flux autocorrelation from heat flux", verbose=verbose)

    kw = {"verbose": False, **kw_correlate}

    flux = array
    ff = flux.dropna(dims.time)
    drift = prefactor * flux.mean(axis=0)
    da = get_autocorrelationNd(ff, off_diagonal=True, **kw)

    if assert_vanishing_mean:
        assert np.mean(drift) < 1e-12, drift

    # prefactor
    da *= prefactor

    da.attrs = {keys.gk_prefactor: prefactor}
    timer()

    return da


def get_heat_flux_power_spectrum(
    flux: xr.DataArray = None, prefactor: float = 1.0, verbose: bool = True
) -> xr.DataArray:
    """compute heat flux power spectrum

    Args:
        flux (xr.DataArray, optional): heat flux [Nt, 3. Defaults to xr.DataArray.
        prefactor (float, optional): prefactor to convert to W/mK/fs
        verbose (bool, optional): be verbose

    Returns:
        xr.DataArray: heat_flux_power_spectrum
    """
    return get_power_spectrum(flux)
