""" generalized Fourier Transforms from positions in supercells to normal modes and
    back """

import numpy as np
from hilde.helpers.lattice_points import map_L_to_i


def u_s_to_u_I(u_q, q_points, lattice_points, eigenvectors, indeces):
    r""" u_iL = 1/sqrt(N) \sum_(q,s) \exp(i q.R_L) e_is(q) u_s(q)

    REM shapes:
        eigenvectors.shape = [n_q, n_s, n_s] """

    n_atoms = len(indeces)
    L_maps = map_L_to_i(indeces)

    u_I = np.zeros([n_atoms, 3])

    for LL, R_L in enumerate(lattice_points):
        # make product of quantities with shape [N_q, N_s, N_s]
        u_temp = (
            np.exp(2j * np.pi * q_points @ R_L)[:, None, None]
            * eigenvectors
            * u_q[:, None, :]
        )
        # sum and reshape to [N_atoms, 3]
        u_I[L_maps[LL]] = (u_temp).sum(axis=(0, 2)).reshape(-1, 3)

    # normalize 1/sqrt(N)
    u_I /= len(q_points) ** 0.5
    return np.array(u_I).real


def u_I_to_u_s(u_I, q_points, lattice_points, eigenvectors, indeces):
    r""" u_s(q) = 1/sqrt(N) e_is(q) . \sum_iL \exp(-i q.R_L) u_iL """

    n_q, n_s = eigenvectors.shape[0:2]
    L_maps = map_L_to_i(indeces)

    u_s = np.zeros([n_q, n_s])

    u_L = np.zeros([n_q, n_s], dtype=complex)

    for LL, R_L in enumerate(lattice_points):
        u_L += (
            np.exp(-2j * np.pi * q_points @ R_L)[:, None]
            * u_I[L_maps[LL]].flatten()[None, :]
        )

    # swapaxes effectively transposes eigenvectors at each q_point
    ievs = eigenvectors.conj().swapaxes(1, 2)

    u_s = (ievs * u_L[:, None, :]).sum(axis=2)

    # normalize 1/sqrt(N)
    u_s /= len(q_points) ** 0.5

    return np.array(u_s).real


def get_A_qst2(in_U_t, in_V_t, in_omegas2):
    r""" compute squared amplitude from mass scaled positions and velocities

        A^2_s(q, t) = u^2_s(q, t) + \omega_s(q)**-2 * \dot{u}^2_s(q, t)

        Parameters:

        in_U_t: list [N_t, N_atoms, 3]
            mass scaled displacements for each time step
        in_V_t: list [N_t, N_atoms, 3]
            mass scaled velocities for each time step
        in_omegas2: list [N_q, N_s]
            eigenvalues (= squared frequencies) of dynamical matrices at commensurate
            q-points
        """

    U_t = np.array(in_U_t)
    V_t = np.array(in_V_t)

    omegas2 = np.array(in_omegas2)
    omegas2[0, :3] = 1e12

    A_qst2 = abs(U_t) ** 2 + omegas2[None, :, :] ** -1 * abs(V_t) ** 2

    A_qst2[:, 0, :3] = 0

    return A_qst2


def get_phi_qst(in_U_t, in_V_t, in_omegas, in_times=None):
    r""" compute phases from mass scaled positions and velocities

    phi_s(q, t) = atan2( -\dot{u}_s(q, t) / \omega_s(q, t) / u_s(q, t))

    Parameters:

    in_U_t: list [N_t, N_atoms, 3]
        mass scaled displacements for each time step
    in_V_t: list [N_t, N_atoms, 3]
        mass scaled velocities for each time step
    in_omegas: list [N_q, N_s]
        frequencies from dynamical matrices at commensurate q-points

    """

    U_t = np.array(in_U_t)
    V_t = np.array(in_V_t)

    omegas = np.array(in_omegas)
    omegas[0, :3] = -1

    if in_times is None:
        omega_t = np.zeros_like(V_t)
    else:
        times = np.array(in_times)
        omega_t = omegas[None, :, :] * times[:, None, None]

    phi_qst = np.arctan2(- V_t - omega_t, omegas[None, :, :] * U_t - omega_t)

    # phase not well defined for 0 modes, set to 0:
    phi_qst[:, 0, :3] = 0

    return phi_qst


def get_E_qst(in_U_t, in_V_t, in_omegas2):
    """ compute mode resolved energies from mass scaled positions and velocities """

    omegas2 = np.array(in_omegas2)

    A_qst2 = get_A_qst2(in_U_t, in_V_t, omegas2)

    E_qst = 0.5 * omegas2[None, :, :] * A_qst2

    return E_qst
