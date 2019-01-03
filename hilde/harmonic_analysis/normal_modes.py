""" generalized Fourier Transforms from positions in supercells to normal modes and
    back """

import numpy as np
from hilde.helpers.lattice_points import map_L_to_i


def u_s_to_u_I(u_q, q_points, lattice_points, eigenvectors, indeces):
    r""" u_iL = 1/N \sum_(q,s) \exp(i q.R_L) e_is(q) u_s(q)

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
    r""" u_s(q) = \sum_iL \exp(-i q.R_L) e_is(q) u_iL """

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
