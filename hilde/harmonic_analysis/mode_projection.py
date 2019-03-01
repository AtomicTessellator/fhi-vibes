""" high-level access to mode projection functionality """

import numpy as np
import scipy.linalg as la

from ase import Atoms

from hilde.helpers.lattice_points import (
    map_I_to_iL,
    get_lattice_points,
    get_commensurate_q_points,
)
from hilde.helpers import Timer
from hilde.helpers.numerics import clean_matrix
from hilde.structure.misc import get_sysname
from hilde.spglib.q_mesh import get_ir_reciprocal_mesh
from .dynamical_matrix import get_dynamical_matrices
from .displacements import get_U, get_dUdt
from .normal_modes import u_I_to_u_s, get_A_qst2, get_phi_qst

# from .force_constants import reshape_force_constants


class HarmonicAnalysis:
    """ provide tools to perform harmonic analysis in periodic systems """

    def __init__(
        self, primitive, supercell, force_constants=None, q_points=None, verbose=False
    ):
        """ Initialize lattice points, commensurate q-points, and solve eigenvalue
            problem at each q-point. Save the results """

        timer = Timer(f"Set up harmonic analysis for {get_sysname(primitive)}:")

        vbsty = {"verbose": verbose}

        self.primitive = primitive
        self.supercell = supercell
        self.force_constants = force_constants

        # intialize
        self._dynamical_matrices = None
        self._irreducible_q_points = None
        self._irreducible_q_points_mapping = None

        # find lattice points:
        self.lattice_points, _ = get_lattice_points(primitive, supercell, **vbsty)

        # find commensurate q_points
        if q_points is None:
            self.q_points = get_commensurate_q_points(primitive, supercell, **vbsty)
        else:
            self.q_points = q_points

        if self.force_constants is None:
            print(f"** Force constants not set, your choice.")
            timer()
            return

        # Write as property instead
        ## solve eigenvalue problem
        # self.omegas2, self.eigenvectors = self.diagonalize_dynamical_matrices()

        ## square root respecting the sign
        # self.omegas = np.sign(self.omegas2) * np.sqrt(abs(self.omegas2))

        # # find map from supercell to primitive + lattice point
        # self.indeces = map_I_to_iL(primitive, supercell)

        timer()

    @property
    def q_points_frac(self):
        """ return relative q points
        Fractional q points: frac_q = q . L.T """

        return np.round(self.q_points @ self.supercell.cell.T).astype(int)

    @property
    def q_points_frac_primitive(self):
        """ return relative q points
        Fractional q points: frac_q = q . L.T """

        return clean_matrix(self.q_points @ self.primitive.cell.T)

    def set_irreducible_q_points(self, is_time_reversal=True, symprec=1e-5):
        """ determine the irreducible q grid in fractionals + mapping """

        mapping, ir_grid = get_ir_reciprocal_mesh(
            self.q_points_frac,
            self.primitive,
            self.supercell,
            is_time_reversal=is_time_reversal,
            symprec=symprec,
        )

        self._irreducible_q_points_mapping = mapping
        self._irreducible_q_points = ir_grid

    def get_irreducible_q_points_frac(self, is_time_reversal=True, symprec=1e-5):
        """ return the irreducible q grid in fractionals + mapping """

        if self._irreducible_q_points is None:
            self.set_irreducible_q_points(
                is_time_reversal=is_time_reversal, symprec=symprec
            )

        return self._irreducible_q_points

    def get_irreducible_q_points_mapping(self, is_time_reversal=True, symprec=1e-5):
        """ return the map from q points to irreducibe qpoints """
        if self._irreducible_q_points is None:
            self.set_irreducible_q_points(
                is_time_reversal=is_time_reversal, symprec=symprec
            )

        return self._irreducible_q_points_mapping

    @property
    def irreducible_q_points_mapping(self):
        """ return mapping from full to irred. q points """
        return self.get_irreducible_q_points_mapping()

    @property
    def irreducible_q_points_frac(self):
        """ return irreducible qpoints in basis of the reciprocal lattice """
        return self.get_irreducible_q_points_frac()

    @property
    def irreducible_q_points(self):
        """ return irreducible qpoints in cartesian """
        ir_grid = self.get_irreducible_q_points_frac()

        return clean_matrix(ir_grid @ self.supercell.get_reciprocal_cell())

    @property
    def irreducible_q_points_frac_primitive(self):
        """ return irreducible qpoints in basis of the reciprocal lattice """
        ir_grid = self.irreducible_q_points

        return clean_matrix(ir_grid @ self.primitive.cell.T)

    def get_dynamical_matrices(self, q_points=None):
        """ return dynamical matrices at (commensurate) q-points """

        if q_points is not None:
            return get_dynamical_matrices(
                q_points, self.primitive, self.supercell, self.force_constants
            )

        if self._dynamical_matrices is None:
            self._dynamical_matrices = get_dynamical_matrices(
                self.q_points, self.primitive, self.supercell, self.force_constants
            )

        return self._dynamical_matrices

    def diagonalize_dynamical_matrices(self, q_points=None):
        """ solve eigenvalue problem for dyn. matrices at (commensurate) q-points """

        omegas2, eigenvectors = [], []
        for dyn_matrix in self.get_dynamical_matrices(q_points):
            w_2, ev = la.eigh(dyn_matrix)
            omegas2.append(w_2)
            eigenvectors.append(ev)

        omegas2 = np.array(omegas2)
        eigenvectors = np.array(eigenvectors)

        return omegas2, eigenvectors

    def omegas2(self, q_points=None):
        """ return square of angular frequencies """
        omegas2, _ = self.diagonalize_dynamical_matrices(q_points)
        return omegas2

    def omegas(self, q_points=None):
        """ return angular frequencies """
        omegas2 = self.omegas2(q_points)

        ## square root respecting the sign
        return np.sign(omegas2) * np.sqrt(abs(omegas2))

    def project(self, trajectory, atoms0=None, times=None):
        """ perform mode projection for atoms objects in trajectory

        Return:
        Amplitdues, Angles, Energies in shape [q, s, t]
        """

        timer = Timer()

        if atoms0 is None:
            atoms0 = self.supercell

        if isinstance(trajectory, Atoms):
            trajectory = [trajectory]

        # sanity check:
        # if la.norm(atoms0.positions - trajectory[0].positions) / len(atoms) > 1

        indeces = map_I_to_iL(self.primitive, atoms0)

        omegas2, eigenvectors = self.diagonalize_dynamical_matrices()
        omegas = self.omegas()

        U_t = [
            u_I_to_u_s(
                get_U(atoms0, atoms),
                self.q_points,
                self.lattice_points,
                eigenvectors,
                indeces,
            )
            for atoms in trajectory
        ]

        V_t = [
            u_I_to_u_s(
                get_dUdt(atoms),
                self.q_points,
                self.lattice_points,
                eigenvectors,
                indeces,
            )
            for atoms in trajectory
        ]

        A_qst2 = get_A_qst2(U_t, V_t, omegas2)
        phi_qst = get_phi_qst(U_t, V_t, omegas, in_times=times)

        E_qst = 0.5 * omegas2[None, :, :] * A_qst2

        timer("project trajectory")

        return A_qst2, phi_qst, E_qst
