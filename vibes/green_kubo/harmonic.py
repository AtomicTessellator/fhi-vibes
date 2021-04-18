"""tools to deal with harmonic flux in GK"""
import collections

import numpy as np
import scipy.optimize as so
import xarray as xr
from ase import units
from vibes import dimensions, keys
from vibes.brillouin import get_symmetrized_array
from vibes.correlation import get_autocorrelationNd
from vibes.dynamical_matrix import DynamicalMatrix
from vibes.helpers import Timer, talk
from vibes.integrate import get_cumtrapz
from vibes.konstanten import to_W_mK

from . import get_gk_prefactor_from_dataset
from .interpolation import get_interpolation_data


_prefix = "gk.harmonic"


def _talk(*args, **kwargs):
    return talk(*args, **kwargs, prefix=_prefix)


def _exp(x, tau, y0=1):
    """compute exponential decay function y = y0 * exp(-x / tau)"""
    return y0 * np.exp(-x / tau)


def get_lifetimes(mode_occupation_acf: xr.DataArray, xarray: bool = True) -> np.ndarray:
    """Compute mode lifetimes from occupation autocorrelation functinon by fitting exp.

    Formulation:
        assume <n_s(t)n_s> ~ g_s(t) = e^{-t / tau_s}

    Args:
        mode_occupation_acf: n_tsq as computed in GK workflow as [Nt, Ns, Nq] array
        xarray: return as DataArray, otherwise as plain numpy

    Returns:
        the mode resolved lifetimes as [Ns, Nq] array
    """
    timer = Timer("Get lifetimes by fitting to exponential", prefix=_prefix)

    tx = mode_occupation_acf.time
    dt = tx[1] - tx[0]

    shape = mode_occupation_acf.shape[1:]
    tau_sq = -np.ones(shape)

    for ns, nq in np.ndindex(shape):

        corr = mode_occupation_acf[:, ns, nq]

        # fit exponential where corr > 1/e
        x_full = np.arange(len(corr))
        first_drop = x_full[corr < 1 / np.e].min()
        if first_drop < 2:
            _talk(f"** acf drops fast for s, q: {ns}, {nq} set tau_sq = np.nan")
            tau_sq[ns, nq] = np.nan
            continue

        x = x_full[:first_drop]
        y = corr[:first_drop]

        (tau, y0), _ = so.curve_fit(_exp, x, y)

        tau_sq[ns, nq] = tau * dt  # back to time axis

    if xarray:
        dims = mode_occupation_acf.dims[1:]
        array = xr.DataArray(tau_sq, dims=dims, name=keys.mode_lifetime)
    else:
        array = tau_sq

    timer()

    return array


def get_n_tsq(
    U_tIa: np.ndarray,
    V_tIa: np.ndarray,
    masses: np.ndarray,
    e_sqI: np.ndarray,
    w_inv_sq: np.ndarray,
    stride: int = 1,
) -> np.ndarray:
    """get mode occupation in shape [Nt, Ns, Nq]

    Formulation:
        n_tsq = |a+_tsq| ** 2 + |a-_tsq| ** 2 = 2 * |a_tsq| ** 2,
            where a_tsq are complex amplitudes as defined, e.g., in Eq. 38.39 or 38.40.

    Args:
        U_tIa:    displacements as [time, atom_index, cartesian_index]
        V_tIa:    velocities    as [time, atom_index, cartesian_index]
        masses:   atomic masses as [atom_index]
        e_sqI: mode eigenvector as [band_index, q_point_index, atom_and_cartesian_index]
        w_inv_sq: inverse frequencies as [band_index, q_point_index]
        stride: use every STRIDE timesteps

    Returns:
        n_tsq: time resolved mode occupations

    """
    Nt = len(U_tIa[::stride])  # no. of timesteps
    # mass-weighted coordinates
    u_tI = np.asarray(masses[None, :, None] ** 0.5 * U_tIa[::stride]).reshape([Nt, -1])
    p_tI = np.asarray(masses[None, :, None] ** 0.5 * V_tIa[::stride]).reshape([Nt, -1])

    # project on modes
    u_tsq = np.moveaxis(e_sqI @ u_tI.T, -1, 0)
    p_tsq = np.moveaxis(e_sqI @ p_tI.T, -1, 0)

    # complex amplitudes and return occupations (= amplitude squared)
    a_tsq = 0.5 * (u_tsq + 1.0j * w_inv_sq[None, :] * p_tsq)

    return 2 * abs(a_tsq) ** 2


def get_J_ha_q(
    n_tsq: np.ndarray, w2_sq: np.ndarray, v_sq: np.ndarray, volume: float
) -> np.ndarray:
    """compute harmonic flux in momentum space (diagonal part)

    Args:
        n_tsq: mode occupations [Nt, Ns, Nq]
        w2_sq: squared frequencies (i.e. eigenvalues of dyn. matrix) [Ns, Nq]
        v_sq: group velocities in Cart. coords [Ns, Nq]
        volume: system volume

    Returns:
        J_ha_q: the harmonic flux per time step

    """
    _ = None
    V = float(volume)
    # make sure we're dealing with ndarrays
    v_sq = np.asarray(v_sq)
    w2_sq = np.asarray(w2_sq)
    n_tsq = np.asarray(n_tsq)

    # J = 1/V * w**2 * n(t) * v
    J_ha_q_sq = 1 / V * (w2_sq[_, :] * n_tsq)[:, :, :, _] * v_sq[_, :]  # [t, s, q, a]

    return J_ha_q_sq.sum(axis=(1, 2))  # -> [t, a]


def get_flux_ha_q_data(
    dataset: xr.Dataset, dmx: DynamicalMatrix, stride: int = 1
) -> collections.namedtuple:
    """return harmonic flux data with flux and mode occupation as xr.DataArrays

    Args:
        dataset: the trajectory dataset
        dmx: the dynamical matrix (object)
        stride: time step length

    Returns:
        (J_ha_q (flux), n_tsq (mode occupation))

    """
    n_tsq = get_n_tsq(
        U_tIa=dataset.displacements,
        V_tIa=dataset.velocities,
        masses=dataset.masses,
        e_sqI=dmx.e_sqI,
        w_inv_sq=dmx.w_inv_sq,
        stride=stride,
    )

    J_ha_q = get_J_ha_q(
        n_tsq=n_tsq,
        w2_sq=dmx.w2_sq,
        v_sq=dmx.v_sqa_cartesian,
        volume=dataset.volume.data.mean(),
    )

    # make dataarrays incl. coords and dims
    ta = (keys.time, dimensions.a)
    tsq = (keys.time, dimensions.s, dimensions.q)
    coords = dataset.coords
    n_tsq = xr.DataArray(n_tsq, coords=coords, dims=tsq, name=keys.mode_occupation)
    J_ha_q = xr.DataArray(J_ha_q, coords=coords, dims=ta, name=keys.hf_ha_q)

    data = {ar.name: ar for ar in (J_ha_q, n_tsq)}

    return collections.namedtuple("flux_ha_q_data", data.keys())(**data)


def get_kappa(
    v_sqa: np.ndarray,
    tau_sq: np.ndarray,
    cv_sq: np.ndarray = None,
    weights: np.ndarray = None,
    scalar: bool = False,
    xarray: bool = True,
) -> xr.DataArray:
    """return BTE-type kappa = c v^2 t for given heat capacity, group vel., and lifetime

    Formula:
        kappa^ab = sum_s c_s t_s v^a_s v^b_s

    Args:
        v_sqa: mode group velocities [s, q, a]
        tau_sq: mode lifetimes
        cv_sq: mode heat capacity per volume =  1/V.kB.T * w_sq ** 4 * <n_sq ** 2>
        weights: weight of q-point in symmetry-reduced grids (tensor will be off!)
        scalar: return scalar, otherwise full tensor
        xarray: return as DataArray instead of ndarray

    Returns:
        thermal condcutivity tensor or scalar
    """
    _ = None

    if weights is None:
        weights = 1.0
    if cv_sq is None:
        cv_sq = 1.0

    # make sure we're dealilng with arrays, cv might be scalar
    cv_sq = np.asarray(cv_sq)
    v_sqa = np.asarray(v_sqa)
    v2_sqa = v_sqa[:, :, :, _] * v_sqa[:, :, _, :]
    kappa_sqab = to_W_mK * np.asarray(weights * cv_sq * tau_sq)[:, :, _, _] * v2_sqa

    kappa_ab = kappa_sqab.sum(axis=(0, 1))

    if xarray:
        array = xr.DataArray(kappa_ab, dims=dimensions.ab, name=keys.mode_lifetime)
    else:
        array = kappa_ab

    if scalar:
        return np.diagonal(array).mean()
    return array


def get_gk_ha_q_data(
    dataset: xr.Dataset,
    dmx: DynamicalMatrix = None,
    interpolate: bool = False,
    stride: int = 1,
) -> collections.namedtuple:
    """return GK dataset entries derived from harmonic model

    Args:
        dataset: the trajectory dataset
        dmx: the dynamical matrix (object), otherwise obtain from dataset
        interpolate: interpolate harmonic flux to dense grid
        stride: time step length

    Returns:
        heat_flux_harmonic_q: J_ha_q = 1/V sum_sq w_sq**2 n_tsq v_sq
        heat_flux_harmonic_q_acf: <J_ha_q (t) J_ha_q>
        heat_flux_harmonic_q_acf_integral: cumulative integral of ACF
        mode_occupation: n_tsq = 2 * |a_tsq|**2
        mode_occupation_acf: respective ACF
        mode_heat_capacity: 1/V.kB.T**2 * w_sq**4 * <n_tsq**2> (intensive/per volume)
    """
    if dmx is None:
        dmx = DynamicalMatrix.from_dataset(dataset)

    J_ha_q, n_tsq = get_flux_ha_q_data(dataset=dataset, dmx=dmx)

    # prefactor
    gk_prefactor = get_gk_prefactor_from_dataset(dataset, verbose=False)

    kw = {"verbose": False}
    jhq = J_ha_q - J_ha_q.mean(axis=0)
    J_ha_q_acf = gk_prefactor * get_autocorrelationNd(jhq, off_diagonal=True, **kw)

    K_ha_q_cum = get_cumtrapz(J_ha_q_acf)

    # occupation autocorrelation normalized to 1
    dn_tsq = n_tsq - n_tsq.mean(axis=0)
    dn_tsq.name = n_tsq.name

    # autocorrelation of occupations
    dn_acf = get_autocorrelationNd(dn_tsq / dn_tsq.std(axis=0), **kw)

    # heat capacity
    volume = dataset.volume.mean().data
    temperature = dataset.temperature.mean().data
    cv_sq = dmx.w2_sq ** 2 * n_tsq.var(axis=0) / units.kB / temperature ** 2 / volume
    cv_sq.name = keys.mode_heat_capacity

    # lifetimes
    tau_sq = get_lifetimes(dn_acf)
    # symmetrized
    map2ir, map2full = dmx.q_grid.map2ir, dmx.q_grid.ir.map2full
    tau_symmetrized_sq = get_symmetrized_array(tau_sq, map2ir=map2ir, map2full=map2full)

    # tau_inv_sq = 1 / tau_sq
    # tau_inv_sq.name = keys.mode_lifetime_inverse

    # thermal concuctivity
    v_sqa = dmx.v_sqa_cartesian
    K_ha_q = get_kappa(v_sqa, tau_sq, cv_sq=cv_sq)
    K_ha_q.name = keys.kappa_ha

    # scalar cv
    k = K_ha_q.data.diagonal().mean()
    cv = k / get_kappa(dmx.v_sqa_cartesian, tau_sq, scalar=True)
    cv = xr.DataArray(cv, name=keys.heat_capacity)

    # symmetrized thermal conductivity
    K_ha_q_symmetrized = get_kappa(v_sqa, tau_sq=tau_symmetrized_sq, cv_sq=cv)
    K_ha_q_symmetrized.name = keys.kappa_ha_symmetrized

    arrays = [J_ha_q, J_ha_q_acf, K_ha_q_cum, K_ha_q, K_ha_q_symmetrized]
    arrays += [n_tsq, dn_acf, cv_sq, cv, tau_sq, tau_symmetrized_sq]

    # add dynamical matrix arrays
    arrays += dmx._get_arrays()

    data = {ar.name: ar for ar in arrays}

    # interpolate
    if interpolate:
        results = get_interpolation_data(dmx, tau_symmetrized_sq, cv)
        data.update(results)

    return collections.namedtuple("gk_ha_q_data", data.keys())(**data)
