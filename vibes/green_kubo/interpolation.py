import numpy as np
import scipy.stats as st
import xarray as xr
from scipy.interpolate import LinearNDInterpolator, griddata
from vibes import dimensions, keys
from vibes.dynamical_matrix import DynamicalMatrix
from vibes.helpers import talk
from vibes.helpers.lattice_points import get_unit_grid_extended


_prefix = "gk.interpolation"


def _talk(*args, **kwargs):
    return talk(*args, **kwargs, prefix=_prefix)


def interpolate_to_gamma(
    q_points: np.ndarray, array_sq: np.ndarray, tol: float = 1e-9
) -> np.ndarray:
    """Get values at Gamma from interpolating surrounding values

    Args:
        q_points: the training grid, FIRST q-point is Gamma, [Nq, 3]
        array_sq: array in shape [Ns, Nq] used for interpolating to Gamma
        tol: tolerance for wrapping

    Returns:
        array_sq: array with interpolated values at Gamma
    """
    # check if first q-point is gamma
    assert np.linalg.norm(q_points[0]) < tol

    # make sure new points are in [-0.5, 0.5)
    train_qs = (q_points[1:] + 0.5 + tol) % 1 - tol - 0.5

    q_gamma = np.zeros(3)

    for ns, l in enumerate(array_sq[:, 1:]):
        interpolator = LinearNDInterpolator(train_qs, l)
        l_s_at_gamma = float(interpolator(q_gamma))
        array_sq[ns, 0] = l_s_at_gamma

    return array_sq


def interpolate_to_grid(
    q_points: np.ndarray,
    train_array_sq: np.ndarray,
    train_points: np.ndarray,
    tol: bool = 1e-9,
) -> np.ndarray:
    """Interpolate values to new q-grid

    Args:
        q_points: new q-points [Nq, 3]
        train_array_sq: the training points in [Ns, Nq]
        train_points: training grid, new point must lie in convex hull
        tol: finite zero

    Returns:
        values from train_array_sq interpolated to q_points
    """

    # make sure new points are in [0, 1]
    new_points = (q_points + tol) % 1 - tol

    grid_interpolators = []
    for l in train_array_sq:
        interpolator = griddata(train_points, l, new_points)
        grid_interpolators.append(interpolator)

    return np.array(grid_interpolators)


def get_interpolation_data(
    dmx: DynamicalMatrix, lifetimes: xr.DataArray, cv: float, nq_max: int = 20
) -> dict:
    """interpolate BTE-type thermal conductivity to dense grid

    Args:
        dmx: dynamical matrix object used for interpolating frequencies etc.
        lifetimes: the lifetimes at commensurate q-points
        cv: heat capacity
        nq_max: interpolate to 3*[nq_max,], use as baseline for extrapolating to infinity

    Returns:
        results dictionary: {
            "interpolation_fit_slope": slope of scaling with 1/nq
            "interpolation_fit_intercept": intercept of linear fit
            "interpolation_fit_stderr": stderr of fit
            "interpolation_correction_factor": K_interpolated = correction_factor * K
            "interpolation_correction_factor_err": error of correction factor
            "interpolation_array": DataArray of interpolated K values
        }
    """
    from .harmonic import get_kappa

    # define scaled lifetimes
    l_sq = dmx.w2_sq * lifetimes

    # get value at gamma from interpolating
    l_sq = interpolate_to_gamma(dmx.q_points, l_sq)

    # create training data on extended unit grid [0, 1]
    train_grid = get_unit_grid_extended(dmx.q_points)
    train_l_sq = l_sq[:, train_grid.map2extended]

    # sanity check that interpolator is idempotent
    kw_train = {
        "train_array_sq": train_l_sq,
        "train_points": train_grid.points_extended,
    }
    assert np.allclose(l_sq, interpolate_to_grid(dmx.q_points, **kw_train))

    # get kappa w/o interpolation (same as K_ha_q_symmetrized)
    kappa_ha = get_kappa(dmx.v_sqa_cartesian, tau_sq=lifetimes, cv_sq=cv)

    # interpolate
    nqs = np.arange(4, nq_max + 1, 2)

    Nq_init = len(dmx.q_points)  # number of commensurate q-points
    kappas = np.zeros((len(nqs), 3, 3))
    for ii, nq in enumerate(nqs):
        mesh = (nq, nq, nq)

        grid, solution = dmx.get_mesh_and_solution(mesh, reduced=False)

        ir_l_int_sq = interpolate_to_grid(q_points=grid.ir.points, **kw_train)

        tau_int_sq = ir_l_int_sq[:, grid.ir.map2full] * solution.w_inv_sq ** 2

        Nq_new = len(grid.points)

        Nq_eff = Nq_new / Nq_init

        k = get_kappa(solution.v_sqa_cartesian, tau_int_sq, cv) / Nq_eff

        kappas[ii] = k
        _talk(f"{nq:3d}, Nq_eff = {Nq_eff:6.2f}, kappa = {np.diagonal(k).mean()}")

    kappas = xr.DataArray(kappas, dims=("nq", *dimensions.a_b), coords={"nq": nqs})

    # interpolate to infinity assuming convergence with 1/nq (Riemann sum)
    ks = np.diagonal(kappas, axis1=1, axis2=2).mean(axis=1)
    m, y0, *_, stderr = st.linregress(nqs ** -1.0, y=ks)

    k_ha = np.diagonal(kappa_ha).mean()
    nq = len(dmx.q_points) ** (1 / 3)

    correction = -m / nq
    correction_factor = 1 + correction / k_ha
    correction_factor_err = stderr / nq / k_ha

    k_ha_int = correction_factor * k_ha

    _talk(f"Initial harmonic kappa value:       {k_ha:.3f} W/mK")
    err_str = f"+/- {stderr / nq:.3f}"
    _talk(f"Fit intercept:                      {y0:.3f} W/mK")
    _talk(f"Fit intercept - initial value:      {y0 - k_ha:.3f} {err_str}  W/mK")
    _talk(f"Interpolated harm. kappa:           {k_ha_int:.3f} {err_str} W/mK")
    _talk(f"Correction:                         {correction:.3f} {err_str} W/mK")
    err_str = f"+/- {correction_factor_err:.3f}"
    _talk(f"Correction factor:                  {correction_factor:.3f} {err_str}")

    dims_w = (dimensions.s, dimensions.q_int)
    dims_q = (dimensions.q_int, dimensions.a)
    results = {
        keys.interpolation_fit_slope: m,
        keys.interpolation_fit_intercept: y0,
        keys.interpolation_fit_stderr: stderr,
        keys.interpolation_correction: correction,
        keys.interpolation_correction_factor: correction_factor,
        keys.interpolation_correction_factor_err: correction_factor_err,
        keys.interpolation_kappa_array: kappas,
        keys.interpolation_q_points: (dims_q, grid.points),
        keys.interpolation_w_sq: (dims_w, solution.w_sq),
        keys.interpolation_tau_sq: (dims_w, tau_int_sq),
    }

    return results
