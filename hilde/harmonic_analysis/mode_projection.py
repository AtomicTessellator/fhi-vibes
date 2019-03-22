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

# from .dynamical_matrix import get_dynamical_matrices
from .displacements import get_U, get_dUdt
from .normal_modes import u_I_to_u_s, get_A_qst2, get_phi_qst

# from .force_constants import reshape_force_constants


class HarmonicAnalysis:
    """ provide tools to perform harmonic analysis in periodic systems """

    def __init__(
        self,
        primitive,
        supercell,
        force_constants=None,
        lattice_points_frac=None,
        lattice_points=None,
        q_points=None,
        verbose=False,
    ):
        """ Initialize lattice points, commensurate q-points, and solve eigenvalue
            problem at each q-point. Save the results """

        timer = Timer(f"Set up harmonic analysis for {get_sysname(primitive)}:")

        vbsty = {"verbose": verbose}

        self.primitive = primitive
        self.supercell = supercell

        # intialize
        self._dynamical_matrices = None
        self._irreducible_q_points = None
        self._irreducible_q_points_mapping = None

        # set or find lattice points, Cartesian and Fractional:
        if lattice_points_frac is not None:
            self.lattice_points_frac = lattice_points_frac
            lps = clean_matrix(lattice_points_frac @ self.primitive.cell)
            self.lattice_points = lps
        elif lattice_points is not None:
            self.lattice_points = lattice_points
            lps = clean_matrix(lattice_points @ la.inv(self.primitive.cell))
            self.lattice_points_frac = lps
        else:
            self.lattice_points, _ = get_lattice_points(primitive, supercell, **vbsty)
            lps = clean_matrix(self.lattice_points @ la.inv(self.primitive.cell))
            self.lattice_points_frac = lps

        # Attach force constants in the shape fc[N_L, N_i, 3, N_j, 3]
        if force_constants is not None:
            self.force_constants = force_constants

        # Sanity check for dimensions
        dim = self.force_constants.shape
        assert dim[0] == self.lattice_points.shape[0], (dim, self.lattice_points.shape)
        assert dim[1] == dim[3] == len(self.primitive), dim
        assert dim[2] == dim[4] == 3, dim

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
        """ return fractional q points w.r.t to primitive cell
        Fractional q points: frac_q = q . L.T """

        return clean_matrix(self.q_points @ self.primitive.cell.T)

    def set_irreducible_q_points(self, is_time_reversal=True, symprec=1e-5):
        """ determine the irreducible q grid in fractionals + mapping """

        fq_supercell = np.round(self.q_points @ self.supercell.cell.T).astype(int)

        mapping, ir_grid = get_ir_reciprocal_mesh(
            fq_supercell,
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
    def irreducible_q_points(self):
        """ return irreducible qpoints in cartesian """
        ir_grid = self.get_irreducible_q_points_frac()

        return clean_matrix(ir_grid @ self.supercell.get_reciprocal_cell())

    @property
    def irreducible_q_points_frac(self):
        """ return irreducible qpoints in basis of the reciprocal lattice """
        ir_grid = self.irreducible_q_points

        return clean_matrix(ir_grid @ self.primitive.cell.T)

    def get_Dx(self):
        """ return mass-scaled hessian """
        na = len(self.primitive)
        fc = self.force_constants.reshape([-1, 3 * na, 3 * na])

        m = self.primitive.get_masses().repeat(3)

        return fc / np.sqrt(m[:, None] * m[None, :])

    def get_Dq(self, q=np.array([0.0, 0.0, 0.0]), fractional=True):
        """ return dynamical matrix at q. """
        if fractional:
            phases = np.exp(2j * np.pi * self.lattice_points_frac @ q)
        else:
            phases = np.exp(2j * np.pi * self.lattice_points @ q)
        Dx = self.get_Dx()

        Dq = (phases[:, None, None] * Dx).sum(axis=0)

        # check for hermiticity
        assert la.norm(Dq - Dq.conj().T) < 1e-12

        return Dq

    def solve_Dq(self, q=np.array([0.0, 0.0, 0.0]), fractional=True):
        """ solve eigenvalue problem for dynamical matrix at q """

        Dq = self.get_Dq(q, fractional=fractional)

        w_2, ev = la.eigh(Dq)

        return w_2, ev

    def diagonalize_dynamical_matrices(self, q_points=None):
        """ solve eigenvalue problem for dyn. matrices at (commensurate) q-points """

        if q_points is None:
            q_points = self.q_points_frac

        omegas2, eigenvectors = [], []
        for q in q_points:
            w_2, ev = self.solve_Dq(q=q)
            omegas2.append(w_2)
            eigenvectors.append(ev)

        omegas2 = np.array(omegas2)
        eigenvectors = np.array(eigenvectors)

        return omegas2, eigenvectors

    def omegas2(self, q_points=None):
        """ return square of angular frequencies """
        omegas2, _ = self.diagonalize_dynamical_matrices(q_points)
        return omegas2

    def eigenvectors(self, q_points=None):
        """ return eigenvectors """
        _, eigenvectors = self.diagonalize_dynamical_matrices(q_points)
        return eigenvectors

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
                self.q_points_frac,
                self.lattice_points_frac,
                eigenvectors,
                indeces,
            )
            for atoms in trajectory
        ]

        V_t = [
            u_I_to_u_s(
                get_dUdt(atoms),
                self.q_points_frac,
                self.lattice_points_frac,
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
