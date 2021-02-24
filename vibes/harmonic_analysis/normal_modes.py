""" generalized Fourier Transforms from positions in supercells to normal modes and
    back """

import numpy as np

from vibes.helpers import warn
from vibes.helpers.lattice_points import map_L_to_i


def get_Zqst(in_Uqst, in_Vqst, in_omegas):
    r""" compute squared amplitude from mass scaled positions and velocities

        Z_s(q, t) = \dot{u}_s(q, t) + \im \omega_s(q) u^2_s(q, t)

        -> E_s(q, t) = |Z_s(q, t)|^2

    Args:
        in_Uqst: mass scaled displacements for each time step
        in_Vqst: mass scaled velocities for each time step
        in_omegas: eigenfrequencies of dynamical matrices at commensurate q-points

    Returns:
        Z_qst: The squared amplitude from mass scaled positions and velocities

    """
    Uqst = np.array(in_Uqst)
    Vqst = np.array(in_Vqst)

    omegas = np.array(in_omegas)
    omegas[0, :3] = 0

    Z_qst = Vqst - 1.0j * omegas[None, :, :] * Uqst

    return Z_qst


def get_A_qst2(in_U_qst, in_V_qst, in_omegas2):
    r""" compute squared amplitude from mass scaled positions and velocities

        A^2_s(q, t) = u^2_s(q, t) + \omega_s(q)**-2 * \dot{u}^2_s(q, t)

    Args:
        in_U_qst: mass scaled displacements for each time step
        in_V_qst: mass scaled velocities for each time step
        in_omegas2: eigenvalues of dynamical matrices at commensurate q-points

    Returns:
        A_qst2: The squared amplitude from mass scaled positions and velocities

    """

    U_qst = np.array(in_U_qst)
    V_qst = np.array(in_V_qst)

    omegas2 = np.array(in_omegas2)
    omegas2[0, :3] = 1e12

    Z_qst = get_Zqst(U_qst, V_qst, omegas2 ** 0.5)

    # A_qst2 = abs(U_qst) ** 2 + omegas2[None, :, :] ** -1 * abs(V_qst) ** 2
    A_qst2 = omegas2[None, :, :] ** -1 * abs(Z_qst) ** 2

    A_qst2[:, 0, :3] = 0

    return A_qst2


def get_phi_qst(in_U_t, in_V_t, in_omegas, in_times=None):
    r""" compute phases from mass scaled positions and velocities

    phi_s(q, t) = atan2( -\dot{u}_s(q, t) / \omega_s(q, t) / u_s(q, t))

    Args:
        in_U_t: [N_t, N_atoms, 3] mass scaled displacements for each time step
        in_V_t: [N_t, N_atoms, 3] mass scaled velocities for each time step
        in_omegas: [N_q, N_s] frequencies from dynamical matrices at commensurate q-points

    Returns:
        The phases from mass scaled positions and velocities

    """
    U_t = np.array(in_U_t).real
    V_t = np.array(in_V_t).real

    warn("Cast mode amplitudes to real values, is this correct?", level=1)

    omegas = np.array(in_omegas)
    omegas[0, :3] = -1

    if in_times is None:
        omega_t = np.zeros_like(V_t)
    else:
        times = np.array(in_times)
        omega_t = omegas[None, :, :] * times[:, None, None]

    phi_qst = np.arctan2(-V_t - omega_t, omegas[None, :, :] * U_t - omega_t)

    # phase not well defined for 0 modes, set to 0:
    phi_qst[:, 0, :3] = 0

    return phi_qst


def get_E_qst(in_U_t, in_V_t, in_omegas2):
    """Compute mode resolved energies from mass scaled positions and velocities

    Args:
        in_U_t: Input mode projected displacements
        in_V_t: Input mode projected velocities
        in_omegas2: Input eigenvalues

    Returns:
        The mode resolved energies from mass scaled positions and velocities

    """
    omegas2 = np.array(in_omegas2)

    A_qst2 = get_A_qst2(in_U_t, in_V_t, omegas2)

    E_qst = 0.5 * omegas2[None, :, :] * A_qst2

    return E_qst


# old stuff, deprecated, but used for testing. why?
def u_I_to_u_s(u_I, q_points, lattice_points, eigenvectors, indeces):
    r""" u_s(q) = 1/sqrt(N) e_is(q) . \sum_iL \exp(-i q.R_L) u_iL

    Args:
        u_I: u in the primitive cell
        q_points: (commensurate) q points
        lattice_points: lattice points within supercell
        eigenvectors: set of eigenvectors of dynamical matrix for each q points
        indeces: index map from supercell index I to image and lattice point index(i, L)

    Returns:
        1/sqrt(N) e_is(q) . \sum_iL \exp(-i q.R_L) u_iL

    """
    n_q, n_s = eigenvectors.shape[0:2]
    L_maps = map_L_to_i(indeces)

    u_s = np.zeros([n_q, n_s])

    u_L = np.zeros([n_q, n_s], dtype=complex)

    for LL, R_L in enumerate(lattice_points):
        u_L += (
            np.exp(-2j * np.pi * q_points @ R_L)[:, None]
            * u_I[L_maps[LL]].flatten()[None, :]
        )

    # assert it's real
    # assert np.linalg.norm(u_L.imag) < 1e-14, (LL, R_L, u_L)

    # swapaxes effectively transposes eigenvectors at each q_point
    ievs = eigenvectors.conj().swapaxes(1, 2)

    u_s = (ievs * u_L[:, None, :]).sum(axis=2)

    # normalize 1/sqrt(N)
    u_s /= len(q_points) ** 0.5

    return np.array(u_s)
