"""compute and analyze heat fluxes"""
import xarray as xr
from ase import units

from vibes import dimensions as dims
from vibes import keys
from vibes.correlation import get_autocorrelationNd
from vibes.fourier import get_fourier_transformed
from vibes.helpers import xtrace
from vibes.integrate import get_cumtrapz

from .utils import Timer, talk


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


def get_kappa_cumulative_dataset(
    dataset: xr.Dataset,
    full: bool = False,
    aux: bool = False,
    delta: str = "auto",
    verbose: bool = True,
    **kw_correlate,
) -> xr.Dataset:
    """compute heat flux autocorrelation and cumulative kappa from heat_flux_dataset

    Args:
        dataset (xr.Dataset): contains heat flux per atom
        full (bool): return correlation function per atom
        aux (bool): add auxiliary heat flux
        delta: compute mode for avalanche time
        kw_correlate (dict): kwargs for `correlate`
    Returns:
        dataset containing
            hfacf, hfacf_scalar, kappa_cumulative, kappa_cumulative_scalar
    """
    pref = get_gk_prefactor_from_dataset(dataset, verbose=verbose)

    kw = {"prefactor": pref, "verbose": verbose, **kw_correlate}

    if full:
        fluxes = dataset[keys.heat_fluxes].copy()
        if aux:
            fluxes += dataset[keys.heat_fluxes_aux]
        J_corr = get_heat_flux_aurocorrelation(fluxes=fluxes, **kw)
    else:
        flux = dataset[keys.heat_flux].copy()
        if aux:
            flux += dataset[keys.heat_flux_aux]
        J_corr = get_heat_flux_aurocorrelation(flux=flux, **kw)

    J_corr.name = keys.hfacf

    kappa = get_cumtrapz(J_corr)
    kappa.name = keys.kappa_cumulative

    # scalars
    Js = xtrace(J_corr) / 3
    Js.name = keys.hfacf_scalar
    Ks = xtrace(kappa) / 3
    Ks.name = keys.kappa_cumulative_scalar

    datasets = (J_corr, Js, kappa, Ks)
    dct = {da.name: da for da in datasets}
    coords = J_corr.coords
    attrs = J_corr.attrs

    return xr.Dataset(dct, coords=coords, attrs=attrs)


def get_heat_flux_aurocorrelation(
    flux: xr.DataArray = None,
    fluxes: xr.DataArray = None,
    prefactor: float = 1.0,
    verbose: bool = True,
    **kw_correlate,
) -> xr.DataArray:
    """compute heat flux autocorrelation function from heat flux Dataset

    Args:
        flux (xr.DataArray, optional): heat flux [Nt, 3. Defaults to xr.DataArray.
        fluxes (xr.DataArray, optional): atomic heat flux [Nt, Na, 3]
        prefactor (float, optional): prefactor to convert to W/mK/fs
        verbose (bool, optional): be verbose
        kw_correlate (dict): kwargs for `correlate`

    Returns:
        xr.DataArray: heat_flux_autocorrelation in W/mK/fs
    """
    timer = Timer("Get heat flux autocorrelation from heat flux", verbose=verbose)

    kw = {"verbose": False, **kw_correlate}

    if fluxes is not None:
        talk(".. return full tensor [N_t, N_a, N_a, 3, 3]", verbose=verbose)
        ff = fluxes.dropna(dims.time)
        da = get_autocorrelationNd(ff, off_diagonal_atoms=True, **kw)
    else:
        ff = flux.dropna(dims.time)
        da = get_autocorrelationNd(ff, off_diagonal_coords=True, **kw)

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
    kw = {"verbose": verbose}
    timer = Timer("Get heat flux power spectrum from heat flux", **kw)

    Jw = get_fourier_transformed(flux.dropna(dims.time), **kw)

    Jspec = abs(Jw)

    timer()

    return Jspec
