"""tools to deal with harmonic flux in GK"""
import collections

import numpy as np
import xarray as xr
from vibes import dimensions, keys
from vibes.correlation import get_autocorrelationNd
from vibes.dynamical_matrix import DynamicalMatrix
from vibes.integrate import get_cumtrapz

from . import get_gk_prefactor_from_dataset


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
        v_sq=dmx.v_sq_cartesian,
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


def get_gk_ha_q_data(
    dataset: xr.Dataset, dmx: DynamicalMatrix = None, stride: int = 1
) -> collections.namedtuple:
    """return GK dataset entries derived from harmonic model

    Args:
        dataset: the trajectory dataset
        dmx: the dynamical matrix (object), otherwise obtain from dataset
        stride: time step length

    Returns:
        (J_ha_q (flux), n_tsq (mode occupation), J_ha_q_cumtrapz (kappa))

    """
    if dmx is None:
        dmx = DynamicalMatrix.from_dataset(dataset)

    J_ha_q, n_tsq = get_flux_ha_q_data(dataset=dataset, dmx=dmx)

    # prefactor
    gk_prefactor = get_gk_prefactor_from_dataset(dataset, verbose=False)

    kw = {"verbose": False}
    jhq = J_ha_q - J_ha_q.mean(axis=0)
    J_ha_q_acf = gk_prefactor * get_autocorrelationNd(jhq, off_diagonal=True, **kw)

    K_ha_q = get_cumtrapz(J_ha_q_acf)

    # occupation autocorrelation normalized to 1
    dn_tsq = n_tsq - n_tsq.mean(axis=0)
    dn_tsq.name = n_tsq.name

    # autocorrelation of occupations
    dn_acf = get_autocorrelationNd(dn_tsq / dn_tsq.std(axis=0), **kw)

    data = {ar.name: ar for ar in (J_ha_q, J_ha_q_acf, K_ha_q, n_tsq, dn_acf)}

    return collections.namedtuple("gk_ha_q_data", data.keys())(**data)
