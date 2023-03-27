"""Green-Kubo post processing"""
from collections import namedtuple

import numpy as np
import xarray as xr
from ase import units
from scipy import signal as sl
from vibes import defaults
from vibes import dimensions as dims
from vibes import keys
from vibes.correlation import get_autocorrelationNd
from vibes.fourier import get_fourier_transformed
from vibes.helpers import Timer, talk, warn
from vibes.helpers.filter import get_filtered
from vibes.integrate import get_cumtrapz
from vibes.konstanten import to_W_mK


_prefix = "GreenKubo"

Timer.prefix = _prefix


def _talk(msg, **kw):
    """wrapper for `utils.talk` with prefix"""
    return talk(msg, prefix=_prefix, **kw)


def gk_prefactor(volume: float, temperature: float, verbose: bool = False) -> float:
    """compute GK prefactor V / (k_B * T^2) and convert from eV/AA^2/fs to W/mK

    Args:
        volume (float): volume of the supercell in AA^3
        temperature (float): avg. temp. in K (trajectory.temperatures.mean())

    Returns:
        V / (k_B * T^2) * 1.602e6
    """
    V = float(volume)
    T = float(temperature)
    prefactor = 1 / units.kB / T ** 2 * V * to_W_mK
    msg = [
        f"Compute Prefactor:",
        f".. Volume:        {V:10.2f}  AA^3",
        f".. Temperature:   {T:10.2f}  K",
        f"-> Prefactor:     {prefactor:10.3e}  W/mK / (eV/AA^2/fs)",
    ]
    _talk(msg, verbose=verbose)
    return float(prefactor)


def get_gk_prefactor_from_dataset(dataset: xr.Dataset, verbose: bool = True) -> float:
    """get the GK prefactor for the dataset, wraps `gk_prefactor`"""
    volume = dataset.attrs[keys.volume]
    temperature = dataset[keys.temperature].mean()
    return gk_prefactor(volume=volume, temperature=temperature, verbose=verbose)


def get_hf_data(
    flux: xr.DataArray,
    dropna_dim=keys.time,
    distribute: bool = True,
    prefactor: float = 1.0,
) -> namedtuple:
    """Compute heat flux autocorrelation and integrated kappa from heat flux

        Args:
            flux [N_t, 3]: the heat flux in an xr.DataArray
            dropna_dim: drop nan values along this dimension (default: `time`)
            distribute: use multiprocessing to parallelize autocorrelation
            prefactor: GK prefactor  V/kB/T**2

    Returns:
        namedtuple: (heat_flux_acf, heat_flux_acf_integral)
    """
    # dropna
    flux = flux.dropna(dropna_dim)

    # compute the time mean flux of flux <J>
    flux_avg = flux.mean(axis=0)

    # compute heat flux autocorrelation function (HFACF) <J(t)J>
    flux_corr = get_autocorrelationNd(
        flux - flux_avg, off_diagonal=True, verbose=False, distribute=distribute
    )

    # use gk prefactor
    flux_corr *= prefactor

    # get integrated kappa
    kappa = get_cumtrapz(flux_corr)

    # return namedtuple w/ HFACF and integrated kappa
    cls = namedtuple("gk_raw_data", (keys.heat_flux_acf, keys.kappa_cumulative))

    return cls(flux_corr, kappa)


def get_lowest_vibrational_frequency(
    velocities: xr.DataArray,
    threshold_freq: float = 0.1,
    prominence: float = 0.2,
    freq_key: str = keys.omega,
    remove_offset: bool = True,
    verbose: bool = False,
) -> float:
    """get the lowest significant vibrational density from VDOS by peak analysis

    Args:
        velocities [N_t, N_a, 3]: DataArray with atomic velocites
        threshold:_freq: neglect data up to this freq in THz (default: 0.1 THz)
        prominence: required prominence for `scipy.signal.find_peaks`
        freq_key: name of the frequency axis (default: `omega`)
        remove_offset: remove offset from VDOS
        verbose: print additional information

    Returns:
        float: the lowest significant vibrational frequency in THz

    """
    # get VDOS (= fourier transform of velocity autocorrelation)
    v = velocities.dropna("time")
    velocties_corr = get_autocorrelationNd(v, verbose=verbose).sum(axis=(1, 2))

    vdos = get_fourier_transformed(velocties_corr, npad=10000, verbose=verbose).real

    if remove_offset:
        vdos -= vdos.min()

    # normalize considering threshold freq.
    vdos /= vdos[vdos[freq_key] > threshold_freq].max()

    # find peaks
    peaks, _ = sl.find_peaks(vdos, prominence=prominence)

    # convert to frequencies
    freqs = vdos[freq_key][peaks]

    # check lowest peak and zero freq.
    v0 = float(vdos[0])
    if v0 > 0.1:
        warn(f"normalized VDOS at omega -> 0 is {v0:.3f}. CHECK?", level=1)

    freq = freqs[0]
    if freq < threshold_freq:
        if freqs[1] > threshold_freq:
            freq = freqs[1]
        else:
            msg = f"Significant peaks are below threshold of {threshold_freq} THZ.\n"
            msg = +f"These Peaks have been found in VDOS (THz), please CHECK!\n,{freqs}"
            warn(msg, level=2)

    return float(freq)


def get_gk_dataset(
    dataset: xr.Dataset,
    interpolate: bool = False,
    window_factor: int = defaults.window_factor,
    filter_prominence: float = defaults.filter_prominence,
    discard: int = 0,
    total: bool = False,
    verbose: bool = True,
) -> xr.Dataset:
    """get Green-Kubo data from trajectory dataset

    Args:
        dataset: a dataset containing `heat_flux` and describing attributes
        interpolate: interpolate harmonic flux to dense grid
        window_factor: factor for filter width estimated from VDOS (default: 1)
        filter_prominence: prominence for peak detection
        discard: discard this many timesteps from the beginning of the trajectory
        total: postprocess gauge-invariant terms of heat flux as well

    Returns:
        xr.Dataset: the processed data

    Workflow:
        1. get heat flux autocorrelation function (HFACF) and the integrated kappa
        2. get lowest significant vibrational frequency and get its time period
        3. filter integrated kappa with this period
        4. get HFACF by time derivative of this filtered kappa
        5. filter the HFACF with the same filter
        6. estimate cutoff time from the decay of the filtered HFACF
        7. run harmonic heat flux and interpolation
    """
    from .harmonic import get_gk_ha_q_data

    # 1. get HFACF and integrated kappa
    heat_flux = dataset[keys.heat_flux]

    if total:  # add non-gauge-invariant contribution
        heat_flux += dataset[keys.heat_flux_aux]

    # get prefactor V/kB/T**2
    gk_prefactor = get_gk_prefactor_from_dataset(dataset, verbose=verbose)

    # get heatflux and integrated kappa
    hfacf, kappa = get_hf_data(heat_flux, prefactor=gk_prefactor)

    # 2. get lowest significant frequency (from VDOS) in THz
    kw = {"prominence": filter_prominence}
    freq = get_lowest_vibrational_frequency(dataset[keys.velocities], **kw)

    # window in fs from freq.:
    window_fs = window_factor / freq * 1000

    kw_talk = {"verbose": verbose}
    _talk(f"Estimate filter window size", **kw_talk)
    _talk(f".. lowest vibrational frequency: {freq:.4f} THz", **kw_talk)
    _talk(f".. corresponding window size:    {window_fs:.4f} fs", **kw_talk)
    _talk(f".. window multiplicator used:    {window_factor:.4f} fs", **kw_talk)

    # 3. filter integrated HFACF (kappa) with this window respecting antisymmetry in time
    kw = {"window_fs": window_fs, "antisymmetric": True, "verbose": verbose}
    k_filtered = get_filtered(kappa, **kw)

    # 4. get the respective HFACF by differentiating w.r.t. time and filtering again
    k_gradient = kappa.copy()

    # compute derivative with j = dk/dt = dk/dn * dn/dt = dk/dn / dt
    dt = float(kappa.time[1] - kappa.time[0])
    k_gradient.data = np.gradient(k_filtered, axis=0) / dt

    j_filtered = get_filtered(k_gradient, window_fs=window_fs, verbose=False)

    # 5. get cutoff times from j_filtered and save respective kappas
    ts = np.zeros([3, 3])
    ks = np.zeros([3, 3])

    for (ii, jj) in np.ndindex(3, 3):
        j = j_filtered[:, ii, jj]
        times = j.time[j < 0]
        if len(times) > 1:
            ta = times.min()
        else:
            warn(f"no cutoff time found", level=1)
            ta = 0
        ks[ii, jj] = k_filtered[:, ii, jj].sel(time=ta)
        ts[ii, jj] = ta

    # report
    if verbose:
        k_diag = np.diag(ks)
        _talk(["Cutoff times (fs):", *np.array2string(ts, precision=3).split("\n")])
        _talk(f"Kappa is:       {np.mean(k_diag):.3f} +/- {np.std(k_diag) / 3**.5:.3f}")
        _talk(["Kappa^ab is: ", *np.array2string(ks, precision=3).split("\n")])

    # 6. compile new dataset
    # add filter parameters to attrs
    attrs = dataset.attrs.copy()
    u = {
        keys.gk_window_fs: window_fs,
        keys.gk_prefactor: gk_prefactor,
        keys.filter_prominence: filter_prominence,
    }
    attrs.update(u)

    data = {
        keys.heat_flux: heat_flux,
        keys.hf_acf: hfacf,
        keys.hf_acf_filtered: j_filtered,
        keys.kappa_cumulative: kappa,
        keys.kappa_cumulative_filtered: k_filtered,
        keys.kappa: (dims.tensor, ks),
        keys.time_cutoff: (dims.tensor, ts),
    }

    # 7. add properties derived from harmonic model
    if keys.fc in dataset:
        data_ha = get_gk_ha_q_data(dataset, interpolate=interpolate)
        data.update(data_ha._asdict())
        data.update({keys.fc: dataset[keys.fc]})

    # add thermodynamic properties
    data.update({key: dataset[key] for key in (keys.volume, keys.temperature)})

    return xr.Dataset(data, coords=kappa.coords, attrs=attrs)
