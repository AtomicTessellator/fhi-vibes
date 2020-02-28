"""compute and analyze heat fluxes"""
import numpy as np
import xarray as xr

# import xrcache as xc
from ase import units
from vibes import defaults
from vibes import dimensions as dims
from vibes import keys
from vibes.correlation import get_autocorrelationNd
from vibes.fourier import get_power_spectrum
from vibes.helpers.xarray import xtrace
from vibes.integrate import get_cumtrapz

from .utils import Timer
from .utils import _talk as talk
from .utils import get_avalanche_data


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
    full: bool = False,
    aux: bool = False,
    delta: str = "auto",
    verbose: bool = True,
    discard: int = 0,
    Fmax: int = defaults.Fmax,
    **kw_correlate,
) -> xr.Dataset:
    """compute heat flux autocorrelation and cumulative kappa from heat_flux_dataset

    Args:
        dataset (xr.Dataset): contains heat flux per atom
        full (bool): return correlation function per atom
        aux (bool): add auxiliary heat flux
        delta: compute mode for avalanche time
        discard: discard this many time steps
        fmax: cutoff for avalanche function
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
    if full:
        talk("** not yet fully implemented, return partially summed heat flux")
        fluxes = dataset[keys.heat_fluxes]
        J_corr = get_jcorr(fluxes).sum(axis=3)
        if aux:
            fluxes = dataset[keys.heat_fluxes] + dataset[keys.heat_fluxes_aux]
            fluxes.name = keys.heat_fluxes_aux
            J_corr_aux = get_jcorr(fluxes)

    else:
        J_corr = get_jcorr(dataset[keys.heat_flux])
        if aux:
            flux = dataset[keys.heat_flux] + dataset[keys.heat_flux_aux]
            flux.name = keys.heat_flux_aux
            J_corr_aux = get_jcorr(flux)

    # compute cumulative kappa
    kappa = get_cumtrapz(J_corr)
    dataarrays = [J_corr, kappa]

    # add aux. heat flux if requested
    if aux:
        dataarrays += [J_corr_aux]

    # create dataset
    dct = {da.name: da for da in dataarrays}
    coords = J_corr.coords
    attrs = J_corr.attrs

    # (scalar) avalanche function
    if full:
        Js = xtrace(J_corr.sum(axis=(1))) / 3
        Ks = xtrace(kappa.sum(axis=(1))) / 3
        kappa = kappa.sum(axis=(1, 3))
    else:
        Js = xtrace(J_corr) / 3
        Ks = xtrace(kappa) / 3

    f_data = get_avalanche_data(Js.to_series(), Fmax=Fmax, verbose=verbose)
    dct.update({keys.avalanche_function: f_data.pop(keys.avalanche_function)})
    attrs.update(f_data)

    # report
    if verbose:
        t_ps = f_data[keys.time_avalanche] / 1000
        n_avalanche = f_data[keys.avalanche_index]
        talk(f"Avalanche time: {t_ps:.3f} ps ({n_avalanche} datapoints)")
        talk(f"Kappa is:       {np.asarray(Ks[n_avalanche])}")
        talk(f"Kappa^ab is: \n{np.asarray(kappa[n_avalanche])}")

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
