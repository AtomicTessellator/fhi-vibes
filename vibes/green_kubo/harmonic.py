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
from vibes.helpers import Timer, progressbar, talk
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


def get_lifetimes(
    mode_energy_acf: xr.DataArray,
    correct_finite_time: bool = True,
    xarray: bool = True,
) -> np.ndarray:
    """Compute mode lifetimes from energy autocorrelation functinon by fitting exp.

    Formulation:
        assume <E_s(t)E_s> / <E^2_s> ~ g_s(t) = e^{-t / tau_s}

    Args:
        mode_energy_acf: E_tsq as computed in GK workflow as [Nt, Ns, Nq] array
        correct_finite_time: add finite-time correction
        xarray: return as DataArray, otherwise as plain numpy

    Returns:
        the mode resolved lifetimes as [Ns, Nq] array
    """
    timer = Timer("Get lifetimes by fitting to exponential", prefix=_prefix)

    tx = mode_energy_acf.time
    dt = tx[1] - tx[0]

    shape = mode_energy_acf.shape[1:]
    tau_sq = -np.ones(shape)

    for ns, nq in np.ndindex(shape):
        corr = mode_energy_acf[:, ns, nq]

        # fit exponential where corr > 1/e
        x_full = np.arange(len(corr))
        x_fit = x_full[corr < 0.5]
        if corr[0] < 1e-12 or x_fit.min() < 2:
            _talk(f"** acf drops fast for s, q: {ns}, {nq} set tau_sq = np.nan")
            tau_sq[ns, nq] = np.nan
            continue

        first_drop = x_fit.min()

        x = x_full[:first_drop]
        y = corr[:first_drop]

        (tau, y0), _ = so.curve_fit(_exp, x, y)

        tau_sq[ns, nq] = tau * dt  # back to time axis

    if correct_finite_time:
        Tmax = mode_energy_acf[keys.time].data.max()
        tau_sq = 1 / (1 / tau_sq - 1 / Tmax)

    if xarray:
        dims = mode_energy_acf.dims[1:]
        array = xr.DataArray(tau_sq, dims=dims, name=keys.mode_lifetime)
    else:
        array = tau_sq

    timer()

    return array


def get_a_tsq(
    U_tIa: np.ndarray,
    V_tIa: np.ndarray,
    masses: np.ndarray,
    e_sqI: np.ndarray,
    w_inv_sq: np.ndarray,
    stride: int = 1,
) -> np.ndarray:
    """get mode amplitude in shape [Nt, Ns, Nq]

    Formulation:
        a_tsq = 1/2 * (u_tsq + 1/iw_sq p_tsq)

    Args:
        U_tIa:    displacements as [time, atom_index, cartesian_index]
        V_tIa:    velocities    as [time, atom_index, cartesian_index]
        masses:   atomic masses as [atom_index]
        e_sqI: mode eigenvector as [band_index, q_point_index, atom_and_cartesian_index]
        w_inv_sq: inverse frequencies as [band_index, q_point_index]
        stride: use every STRIDE timesteps

    Returns:
        a_tsq: time resolved mode amplitudes

    """
    Nt = len(U_tIa[::stride])  # no. of timesteps
    # mass-weighted coordinates
    u_tI = np.asarray(masses[None, :, None] ** 0.5 * U_tIa[::stride]).reshape([Nt, -1])
    p_tI = np.asarray(masses[None, :, None] ** 0.5 * V_tIa[::stride]).reshape([Nt, -1])

    # project on modes
    u_tsq = np.moveaxis(e_sqI @ u_tI.T, -1, 0)
    p_tsq = np.moveaxis(e_sqI @ p_tI.T, -1, 0)

    # complex amplitudes
    a_tsq = 0.5 * (u_tsq + 1.0j * w_inv_sq[None, :] * p_tsq)

    return a_tsq


def get_J_ha_q(E_tsq: np.ndarray, v_sqa: np.ndarray, volume: float) -> np.ndarray:
    """compute harmonic flux in momentum space (diagonal part)

    Args:
        E_tsq: mode-resolved energy [Nt, Ns, Nq]
        w2_sq: squared frequencies (i.e. eigenvalues of dyn. matrix) [Ns, Nq]
        v_sqa: group velocities in Cart. coords [Ns, Nq, 3]
        volume: system volume

    Returns:
        J_ha_q: the harmonic flux per time step

    """
    _ = None
    V = float(volume)
    # make sure we're dealing with ndarrays
    v_sqa = np.asarray(v_sqa)
    E_tsq = np.asarray(E_tsq)

    # J = 1/V * w**2 * n(t) * v
    J_ha_q_sq = 1 / V * (E_tsq)[:, :, :, _] * v_sqa[_, :]  # [t, s, q, a]

    return J_ha_q_sq, J_ha_q_sq.sum(axis=(1, 2))  # -> [t, a]


def get_flux_ha_q_data(
    dataset: xr.Dataset, dmx: DynamicalMatrix, stride: int = 1
) -> collections.namedtuple:
    """return harmonic flux data with flux and mode energies as xr.DataArrays

    Args:
        dataset: the trajectory dataset
        dmx: the dynamical matrix (object)
        stride: time step length

    Returns:
        (J_ha_q (flux), E_tsq (mode energy))

    """
    a_tsq = get_a_tsq(
        U_tIa=dataset.displacements,
        V_tIa=dataset.velocities,
        masses=dataset.masses,
        e_sqI=dmx.e_sqI,
        w_inv_sq=dmx.w_inv_sq,
        stride=stride,
    )

    # compute mode-resolved energy
    # E = 2 * w2 * a+.a
    E_tsq = 2 * abs(a_tsq) ** 2 * dmx.w2_sq
    v_sqa = dmx.v_sqa_cartesian
    volume = dataset.volume.data.mean()

    J_ha_q_sq, J_ha_q = get_J_ha_q(E_tsq=E_tsq, v_sqa=v_sqa, volume=volume)

    # make dataarrays incl. coords and dims
    ta = (keys.time, dimensions.a)
    tsq = (keys.time, dimensions.s, dimensions.q)
    tsqa = (keys.time, dimensions.s, dimensions.q, dimensions.a)
    coords = dataset.coords
    E_tsq = xr.DataArray(E_tsq, coords=coords, dims=tsq, name=keys.mode_energy)
    J_ha_q = xr.DataArray(J_ha_q, coords=coords, dims=ta, name=keys.hf_ha_q)
    J_ha_q_sq = xr.DataArray(J_ha_q_sq, coords=coords, dims=tsqa, name="heat_flux_harmonic_q_tsqa") # tsqa

    data = {ar.name: ar for ar in (J_ha_q_sq, J_ha_q, E_tsq)}

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
    tau_sq = np.nan_to_num(tau_sq)  # replace nan values by 0
    v2_sqa = v_sqa[:, :, :, _] * v_sqa[:, :, _, :]
    kappa_sqab = to_W_mK * np.asarray(weights * cv_sq * tau_sq)[:, :, _, _] * v2_sqa

    kappa_ab = kappa_sqab.sum(axis=(0, 1))

    if xarray:
        array = xr.DataArray(kappa_ab, dims=dimensions.a_b, name=keys.mode_lifetime)
    else:
        array = kappa_ab

    if scalar:
        return np.diagonal(array).mean()
    return array

def full_to_voigt(corr): 
    # map 3x3 tensor to voigt index
    corr += corr.swapaxes(-1,-2)
    corr /= 2
    corr = corr.reshape([*corr.shape[0:-2], 9])
    corr = corr[...,[0,4,8,5,2,1]]
    return corr

def get_mode_resolved_ha_q_acf(dataset, J_ha_q_sq, gk_prefactor, dump_q=True):
    # mode-resolved correlation
    from itertools import product
    from vibes.correlation import c1, get_unique_dimensions
    
    jhq_sq = J_ha_q_sq - J_ha_q_sq.mean(axis=0) # t, s, q, a
    #J_ha_q_sq_acf = gk_prefactor * get_autocorrelationNd(jhq_sq, off_diagonal=True, **kw)
    dims = jhq_sq.dims
    # new_dims = t, s, q, 2, 6 
    _shape = jhq_sq.shape
    new_shape = (_shape[0], 2, 6)
    corr = np.zeros(new_shape)

    jhq_q = jhq_sq.sum(axis=1)
    for _q in range(_shape[2]):
        j2 = jhq_q[:,_q,:] # shape = t, 3
        Nt, *shape = np.shape(j2)
        _j2 = np.asarray(j2).reshape(Nt,-1).swapaxes(0, 1)
        
        q_corr = np.zeros((*_shape[0:-2], 2, 6))
        for _s in range(_shape[1]):
            j1 = jhq_sq[:,_s,_q,:] # shape = t, 3
                
            # reshape to (X, Nt)
            _j1 = np.asarray(j1).reshape(Nt,-1).swapaxes(0, 1)
            
            # compute correlation function
            off_corr = np.zeros([_j1.shape[0] ** 2, Nt])
            for ii, X in enumerate([*product(_j1,_j2)]):
                off_corr[ii] = c1(X)
            off_corr = np.moveaxis(off_corr.reshape((*shape, *shape, Nt)), -1, 0)

            diag_corr = np.zeros([_j1.shape[0] ** 2, Nt])
            for ii, X in enumerate([*product(_j1,_j1)]):
                diag_corr[ii] = c1(X)
            diag_corr = np.moveaxis(diag_corr.reshape((*shape, *shape, Nt)), -1, 0)
            
            off_corr = full_to_voigt(off_corr)
            diag_corr = full_to_voigt(diag_corr)

            q_corr[:,_s,0,:] = gk_prefactor * diag_corr
            q_corr[:,_s,1,:] = gk_prefactor * off_corr

        if dump_q:
            J_ha_q_sq_acf_q = xr.DataArray(
                q_corr,
                dims= (*dims[0:-2], "ifdiag", "voigt"),
                coords=dataset[keys.heat_flux].coords,
                name=f"{jhq_sq.name}_{keys.autocorrelation}",
            )
            J_ha_q_sq_acf_q.astype("f4") # FIXME
            _talk(f"this shape {dataset[keys.heat_flux].coords}")
            _talk(f"our shape now {J_ha_q_sq_acf_q.shape}")
            jha_q_filename = f"J_ha_q_acf_{_q}.nc"
            _talk(f"harmonic heat flux dump to {jha_q_filename}")
            encoding = {J_ha_q_sq_acf_q.name:{"dtype": "f4"}} # FIXME
            J_ha_q_sq_acf_q.to_netcdf(jha_q_filename, encoding=encoding)
        
        corr += q_corr.sum(axis=1)
    _talk(f"corr_shape: {corr.shape}")
    new_dims = (dims[0], "ifdiag", "voigt")
    _talk(f"new_dims: {new_dims}")
    J_ha_q_sq_acf  = xr.DataArray(
        corr,
        dims=new_dims,
        coords=jhq_sq.coords,
        name=f"{jhq_sq.name}_{keys.autocorrelation}",
    )
    return J_ha_q_sq_acf

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
        mode_energy: E_tsq = 2 * w_sq **2 * |a_tsq|**2
        mode_energy_acf: respective ACF
        mode_heat_capacity: 1/V.kB.T**2 * w_sq**4 * <n_tsq**2> (intensive/per volume)
    """
    if dmx is None:
        dmx = DynamicalMatrix.from_dataset(dataset)

    J_ha_q_sq, J_ha_q, E_tsq = get_flux_ha_q_data(dataset=dataset, dmx=dmx)

    # gk prefactor
    gk_prefactor = get_gk_prefactor_from_dataset(dataset, verbose=False)

    kw = {"verbose": False}
    jhq = J_ha_q - J_ha_q.mean(axis=0)
    J_ha_q_acf = gk_prefactor * get_autocorrelationNd(jhq, off_diagonal=True, **kw)

    K_ha_q_cum = get_cumtrapz(J_ha_q_acf)

    # energy fluctuations
    dE_tsq = E_tsq - E_tsq.mean(axis=0)
    dE_tsq.name = E_tsq.name

    # normalize dE_tsq by std. dev. and replace nan values by 0
    dE_norm = dE_tsq / dE_tsq.std(axis=0)
    dE_norm.data = np.nan_to_num(dE_norm)

    # energy autocorrelation
    dE_acf = get_autocorrelationNd(dE_norm, **kw)

    # heat capacity
    volume = dataset.volume.mean().data
    temperature = dataset.temperature.mean().data
    cv_sq = E_tsq.var(axis=0) / units.kB / temperature**2 / volume
    cv_sq.name = keys.mode_heat_capacity

    # get lifetimes from fitting exp. function to mode energy
    tau_sq = get_lifetimes(dE_acf)

    # symmetrize by averaging over symmetry-related q-points
    map2ir, map2full = dmx.q_grid.map2ir, dmx.q_grid.ir.map2full
    tau_symmetrized_sq = get_symmetrized_array(tau_sq, map2ir=map2ir, map2full=map2full)

    # tau_inv_sq = 1 / tau_sq
    # tau_inv_sq.name = keys.mode_lifetime_inverse

    # compute thermal concuctivity
    v_sqa = dmx.v_sqa_cartesian
    K_ha_q = get_kappa(v_sqa=v_sqa, tau_sq=tau_sq, cv_sq=cv_sq)
    K_ha_q.name = keys.kappa_ha

    # scalar cv: account for cv/kB not exactly 1 in the numeric simulation
    # choose such that c * K(c_s=kB) = K(c_s=c_s)
    k = K_ha_q.data.diagonal().mean()
    if k < 1.0e-4:  # comment: this criteria needs a careful tuning
        cv = cv_sq.mean()
    else:
        cv = k / get_kappa(v_sqa=v_sqa, tau_sq=tau_sq, scalar=True)
    cv = xr.DataArray(cv, name=keys.heat_capacity)

    # symmetrized thermal conductivity
    K_ha_q_symmetrized = get_kappa(v_sqa=v_sqa, tau_sq=tau_symmetrized_sq, cv_sq=cv)
    K_ha_q_symmetrized.name = keys.kappa_ha_symmetrized

    arrays = [J_ha_q, J_ha_q_acf, K_ha_q_cum, K_ha_q, K_ha_q_symmetrized]
    arrays += [E_tsq, dE_acf, cv_sq, cv, tau_sq, tau_symmetrized_sq]

    # mode-resolved correlation 
    #timer = Timer("Calc mode-resolved correlation", prefix=_prefix)
    #J_ha_q_sq_acf = get_mode_resolved_ha_q_acf(dataset, J_ha_q_sq, gk_prefactor)
    #K_ha_q_sq_cum = get_cumtrapz(J_ha_q_sq_acf)
    #arrays += [J_ha_q_sq, J_ha_q_sq_acf, K_ha_q_sq_cum]
    #timer()

    # add dynamical matrix arrays
    arrays += dmx._get_arrays()

    data = {ar.name: ar.astype("f4") for ar in arrays} # FIXME

    # interpolate
    if interpolate:
        if keys.interpolation_fc in dataset:
            # use a different FC for interpolation
            from vibes.dynamical_matrix import InterpolationDynamicalMatrix

            _talk("Use another FC for better interpolation!")
            dmx = InterpolationDynamicalMatrix.from_dataset(dataset)
            data.update({keys.interpolation_fc: dataset[keys.interpolation_fc]})

        results = get_interpolation_data(dmx, tau_symmetrized_sq, cv)
        data.update(results)

    return collections.namedtuple("gk_ha_q_data", data.keys())(**data)
