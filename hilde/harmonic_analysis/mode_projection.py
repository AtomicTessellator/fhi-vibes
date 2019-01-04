""" high-level access to mode projection functionality """

import numpy as np
import scipy.linalg as la

from hilde.helpers.lattice_points import (
    map_I_to_iL,
    get_lattice_points,
    get_commensurate_q_points,
)
from hilde.helpers import Timer
from .dynamical_matrix import get_dynamical_matrices
from .displacements import get_U, get_dUdt
from .normal_modes import u_I_to_u_s, get_A_qst2, get_phi_qst

# from .force_constants import reshape_force_constants


class HarmonicAnalysis:
    """ provide tools to perform harmonic analysis in periodic systems """

    def __init__(self, primitive, supercell, force_constants, q_points=None):
        """ Initialize lattice points, commensurate q-points, and solve eigenvalue
            problem at each q-point. Save the results """

        timer = Timer()

        self.primitive = primitive
        self.supercell = supercell
        self.force_constants = force_constants

        # intialize
        self._dynamical_matrices = None

        # find lattice points:
        self.lattice_points, _ = get_lattice_points(primitive, supercell)

        # find commensurate q_points
        if q_points is None:
            self.q_points = get_commensurate_q_points(primitive, supercell)
        else:
            self.q_points = q_points

        # solve eigenvalue problem
        self.omegas2, self.eigenvectors = self.diagonalize_dynamical_matrices()

        # square root respecting the sign
        self.omegas = np.sign(self.omegas2) * np.sqrt(abs(self.omegas2))

        # # find map from supercell to primitive + lattice point
        # self.indeces = map_I_to_iL(primitive, supercell)

        timer(f"setup for harmonic analysis finished")

    @property
    def dynamical_matrices(self):
        """ return dynamical matrices at (commensurate) q-points """

        if self._dynamical_matrices is None:
            self._dynamical_matrices = get_dynamical_matrices(
                self.q_points, self.primitive, self.supercell, self.force_constants
            )

        return self._dynamical_matrices

    def diagonalize_dynamical_matrices(self, q_points=None):
        """ solve eigenvalue problem for dyn. matrices at (commensurate) q-points """

        if q_points is not None:
            self.q_points = q_points

        omegas2, eigenvectors = [], []
        for dyn_matrix in self.dynamical_matrices:
            w_2, ev = la.eigh(dyn_matrix)
            omegas2.append(w_2)
            eigenvectors.append(ev)

        omegas2 = np.array(omegas2)
        eigenvectors = np.array(eigenvectors)

        return omegas2, eigenvectors

    def project(self, trajectory, atoms0=None, times=None):
        """ perform mode projection for atoms objects in trajectory """

        if atoms0 is None:
            atoms0 = self.supercell

        # sanity check:
        # if la.norm(atoms0.positions - trajectory[0].positions) / len(atoms) > 1

        indeces = map_I_to_iL(self.primitive, atoms0)

        U_t = [
            u_I_to_u_s(
                get_U(atoms0, atoms),
                self.q_points,
                self.lattice_points,
                self.eigenvectors,
                indeces,
            )
            for atoms in trajectory
        ]

        V_t = [
            u_I_to_u_s(
                get_dUdt(atoms),
                self.q_points,
                self.lattice_points,
                self.eigenvectors,
                indeces,
            )
            for atoms in trajectory
        ]

        A_qst2 = get_A_qst2(U_t, V_t, self.omegas2)
        phi_qst = get_phi_qst(U_t, V_t, self.omegas, in_times=times)

        E_qst = 0.5 * self.omegas2[None, :, :] * A_qst2

        return A_qst2, phi_qst, E_qst
