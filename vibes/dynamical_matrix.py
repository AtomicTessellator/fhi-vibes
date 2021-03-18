""" Tools for dealing with force constants """
import collections
import numpy as np
from ase import Atoms
from phonopy import Phonopy

from vibes.brillouin import get_q_grid
from vibes.helpers import talk
from vibes.helpers.lattice_points import get_commensurate_q_points
from vibes.io import get_identifier
from vibes.konstanten import gv_to_AA_fs

from .force_constants import ForceConstants
from .spglib import get_ir_reciprocal_mesh


la = np.linalg


_prefix = "dynamical_matrix"


def _talk(msg, verbose=True):
    return talk(msg, prefix=_prefix, verbose=verbose)


solution_keys = ("w_sq", "w_inv_sq", "w2_sq", "v_sq_cartesian", "e_isq")

Solution = collections.namedtuple("solution", solution_keys)


def get_full_solution_from_ir(q_grid, ir_solution):
    """use symmetry to reconstruct the solution at each q from the information at ir. qs

    Args:
        q_grid: grid including the symmetry information
        ir_solution: solution for the irred. grid poins

    Returns:
        solution at each q in the grid

    """
    Ns, Nq = len(ir_solution.w_sq), len(q_grid.points)
    w_sq = np.zeros([Ns, Nq])
    w_inv_sq = np.zeros([Ns, Nq])
    v_sq_cartesian = np.zeros([Ns, Nq, 3])
    w2_sq = np.zeros([Ns, Nq])
    e_isq = np.zeros([Ns, Ns, Nq], dtype=complex)

    for iq, q in enumerate(q_grid.points_cartesian):

        # index of prototype in the full grid
        ip = q_grid.map2ir_points[iq]

        # index of the symmetry op from prototype to q
        ir = q_grid.symop2ir[iq]

        index_map = q_grid.spg_data.index_maps[ir]

        # cartesian rotation matrix
        rot = q_grid.spg_data.rotations_cartesian[ir]

        # check that q = rot . p
        p = q_grid.ir.points_cartesian[ip]
        assert np.allclose(q @ rot, p), (iq, ip, ir, p, q @ rot, rot)

        # symmetry operator matrix representation
        permutation = np.eye(len(index_map))[index_map]
        G = np.kron(permutation, rot)

        # rotate velocities

        # check group velocities
        # vq = rot @ vp = vp @ rot.T
        vp = ir_solution.v_sq_cartesian[:, ip]
        vq = vp @ rot.T

        # rotate eigenvector
        e_isq_rot = G @ ir_solution.e_isq[:, :, ip]

        # assign rotated quantities
        v_sq_cartesian[:, iq] = vq
        e_isq[:, :, iq] = e_isq_rot

        # assign unchanged properties
        w_sq[:, iq] = ir_solution.w_sq[:, ip]
        w_inv_sq[:, iq] = ir_solution.w_inv_sq[:, ip]

    # restore squared frequencies
    w2_sq = np.sign(w_sq) * w_sq ** 2

    return Solution(w_sq, w_inv_sq, w2_sq, v_sq_cartesian, e_isq)


class DynamicalMatrix(ForceConstants):
    """Mass-weighted force constants"""

    def __init__(
        self,
        force_constants: np.ndarray,
        primitive: Atoms,
        supercell: Atoms,
        symmetry: bool = True,
        tol: float = 1e-12,
    ):
        """like ForceConstants, but with mass weighting and q-space things

        Args:
            force_constants: the force constant matrix in phonopy shape ([Np, Ns, 3, 3])
            primitive: the reference primitive structure
            supercell: the reference supercell
            symmetry: use symmetry to reduce q-points

        """
        # init ForceConstants
        super().__init__(
            force_constants=force_constants, primitive=primitive, supercell=supercell,
        )
        # now scale the forceconstants array
        pmasses = self.primitive.get_masses()
        smasses = self.supercell.get_masses()
        self._mass_weights = (pmasses[:, None] * smasses[None, :]) ** -0.5

        # attach phonopy object for computing eigensolutions fast
        smatrix = np.rint(supercell.cell @ primitive.cell.reciprocal().T).T
        phonon = Phonopy(unitcell=primitive, supercell_matrix=smatrix)
        phonon.force_constants = self._fc_phonopy
        self._phonon = phonon

        # attach commensurate q-points (= inverse lattice points)
        pcell, scell = self.primitive.cell, self.supercell.cell
        q_points_commensurate = get_commensurate_q_points(pcell, scell, fractional=True)

        # create a grid including symmetry information
        self._q_grid = get_q_grid(q_points_commensurate, primitive=primitive)

        # solve on the irreducible grid
        kw = {"eigenvectors": True}
        if symmetry:
            self._ir_solution = self.get_solution(q_points=self.q_grid.ir.points, **kw)
            solution = get_full_solution_from_ir(self.q_grid, self.ir_solution)
        else:
            self._ir_solution = None
            solution = self.get_solution(q_points=self.q_grid.points, **kw)
        self._solution = solution

        # map eigenvectors to supercell
        # e_Isq = 1/N**.5 * e^iq.R_I e_isq in vectorized form
        # I: (I, alpha) = supercell atom label and cart. coord
        Rs = self.phonon.supercell.positions
        qs = self.q_points_cartesian
        Nq = len(qs)

        p_indices = self.I2iL_map[:, 0]
        # index map including Cart. coords because i=(i, a), I=(I, a)
        i2I = np.concatenate([np.arange(3) + 3 * ii for ii in p_indices])
        e_isq = self.e_isq[i2I]

        phases = np.exp(2j * np.pi * (qs[None, :] * Rs[:, None, :]).sum(axis=-1))

        self._e_Isq = Nq ** -0.5 * phases.repeat(3, axis=0)[:, None, :] * e_isq

        # sanity check:
        # reconstruct dynamical matrix for supercell from eigensolution
        D_IJ = (
            self.e_Isq[:, None, :, :]
            * self.w2_sq[None, None, :, :]
            * self.e_Isq.conj()[None, :, :, :]
        ).sum(axis=(-1, -2))

        if not np.allclose(D_IJ, self.remapped):
            diff = la.norm(D_IJ - self.remapped)
            talk(f"** dyn. matrix reconstructed from eigensolution differs by {diff}")

        _talk(f"Setup complete, eigensolution is unitary.")

    def __repr__(self):
        dct = get_identifier(self.primitive)
        sg, name = dct["space_group"], dct["material"]
        rep = f"Dynamical Matrix for {name} (sg {sg}) of shape {self.array.shape}"
        return repr(rep)

    @property
    def solution(self):
        """return commensurate solution"""
        return self._solution

    @property
    def ir_solution(self):
        """return commensurate solution on reduced grid"""
        return self._ir_solution

    @property
    def phonon(self) -> Phonopy:
        return self._phonon

    def get_solution(
        self,
        q_points: np.ndarray,
        # weights: np.ndarray = None,
        eigenvectors: bool = False,
    ) -> collections.namedtuple:
        """solve at each q-point, optionally with eigenvectors

        Args:
            q_points: the q points to compute solutions for
            weights: their corresponding weights (if working with symmetry-reduced grid)
            eigenvectors: include eigenvectors

        Return:
            solution
                .w_sq: frequencies
                .w_inv_sq: inverse frequencies respecting zeros
                .w2_sq: squared frequencies (eigenvalues)
                .v_sq_cartesian: group velocities in Cart. coords in (AA/fs)
                # .weights: weight of each q-point
                .e_isq: eigenvectors

        """
        kw = {"with_group_velocities": True, "with_eigenvectors": eigenvectors}
        self.phonon.run_qpoints(q_points, **kw)
        data = self.phonon.get_qpoints_dict()

        # create solution tuple
        w_sq = data["frequencies"].swapaxes(0, 1) / self.phonon._factor
        v_sq_cart = data["group_velocities"].swapaxes(0, 1) / self.phonon._factor
        v_sq_cart *= gv_to_AA_fs
        if eigenvectors:
            e_isq = np.moveaxis(data["eigenvectors"], 0, -1)
        else:
            e_isq = None

        # make 0s 0, take 4th smallest absolute eigenvalue as reference
        thresh = sorted(abs(w_sq.flatten()))[3] * 0.1
        w_sq[abs(w_sq) < thresh] = 0
        w2_sq = np.sign(w_sq) * w_sq ** 2  # eigenvalues

        w = w_sq.copy()
        thresh = sorted(abs(w.flatten()))[3] * 0.9
        w[abs(w) < thresh] = 1e20
        w_inv_sq = 1 / w
        w_inv_sq[abs(w_inv_sq) < 1e-9] = 0

        return Solution(w_sq, w_inv_sq, w2_sq, v_sq_cart, e_isq)

    def get_mesh_and_solution(self, mesh: list) -> tuple:
        """create a qpoints mesh and solutions on the mesh

        Args:
            mesh: mesh numbers, e.g. [4, 4, 4]
        Returns:
            q_grid: the q-points and their weights
            solution: the solution on the mesh

        """
        assert len(mesh) == 3, len(mesh)
        kw = {"atoms": self.primitive, "monkhorst": True}
        ir_mesh = get_ir_reciprocal_mesh(mesh=mesh, **kw)
        ir_solution = self.get_solution(q_points=ir_mesh.points)

        return ir_mesh, ir_solution

    @property
    def q_grid(self):
        """return the grid object representing the commensurate q-point grid"""
        return self._q_grid

    @property
    def q_points(self):
        """return commensurate set of q_points in fractional coords"""
        return self.q_grid.points

    @property
    def q_points_cartesian(self):
        """return commensurate set of q_points in Cartesian coords"""
        return self.q_grid.points_cartesian

    @property
    def v_sq_cartesian(self):
        """cartesian group velocities in [3 * Np, Nq, 3]"""
        return self.solution.v_sq_cartesian

    @property
    def w_sq(self):
        """eigenfrequencies in ASE units [3 * Np, Nq]"""
        return self.solution.w_sq.copy()

    @property
    def w2_sq(self):
        """eigenvalues [3 * Np, Nq]"""
        return self.solution.w2_sq.copy()

    @property
    def w_inv_sq(self):
        """inverse eigenfrequencies [3 * Np, Nq] respecting sign and zeros"""
        return self.solution.w_inv_sq.copy()

    @property
    def e_isq(self):
        """|i, sq> [3 * Np, 3 * Np, Nq]"""
        return self.solution.e_isq.copy()

    @property
    def e_Isq(self):
        """eigenvectors w.r.t. comm. q-points |I, sq> [3 * Na, 3 * Np, Nq]"""
        return self._e_Isq.copy()

    @property
    def e_sqI(self):
        """project on eigenmodes <sq, I| [3 * Np, Nq, 3 * Na]"""
        return np.moveaxis(self.e_Isq, 0, -1)

    @property
    def qij(self):
        """get dynamical matrices at commensurate q-points"""
        D_ijq = (
            self.e_isq[:, None, :, :]
            * self.w2_sq[None, None, :, :]
            * self.e_isq.conj()[None, :, :, :]
        ).sum(axis=-2)

        return np.moveaxis(D_ijq, -1, 0)

    def copy(self):
        """return copy"""
        return DynamicalMatrix(
            force_constants=self._fc_phonopy,
            primitive=self.primitive,
            supercell=self.supercell,
        )
