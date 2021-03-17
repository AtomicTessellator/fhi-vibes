""" Tools for dealing with force constants """

import numpy as np
from ase import Atoms
from phonopy import Phonopy

from vibes.brillouin import get_ir_grid
from vibes.helpers import talk
from vibes.helpers.lattice_points import get_commensurate_q_points
from vibes.io import get_identifier

from .force_constants import ForceConstants


la = np.linalg


_prefix = "dynamical_matrix"


def _talk(msg, verbose=True):
    return talk(msg, prefix=_prefix, verbose=verbose)


class DynamicalMatrix(ForceConstants):
    """Mass-weighted force constants"""

    def __init__(
        self,
        force_constants: np.ndarray,
        primitive: Atoms,
        supercell: Atoms,
        tol: float = 1e-12,
    ):
        """like ForceConstants, but with mass weighting and q-space things

        Args:
            force_constants: the force constant matrix in phonopy shape ([Np, Ns, 3, 3])
            primitive: the reference primitive structure
            supercell: the reference supercell

        """
        # init ForceConstants
        super().__init__(
            force_constants=force_constants, primitive=primitive, supercell=supercell,
        )
        # now scale the forceconstants array
        pmasses = self.primitive.get_masses()
        smasses = self.supercell.get_masses()
        self._mass_weights = (pmasses[:, None] * smasses[None, :]) ** -0.5

        # attach phonopy object
        smatrix = np.rint(supercell.cell @ primitive.cell.reciprocal().T).T
        phonon = Phonopy(unitcell=primitive, supercell_matrix=smatrix)
        phonon.force_constants = self._fc_phonopy
        self._phonon = phonon

        # attach commensurate q-points (= inverse lattice points)
        pcell, scell = self.primitive.cell, self.supercell.cell
        self._commensurate_q_points = get_commensurate_q_points(
            pcell, scell, fractional=True
        )

        self.set_at_commensurate_qs()

        # map eigenvectors to supercell
        # e_Isq = 1/N**.5 * e^iq.R_I e_isq in vectorized form
        # I: (I, alpha) = supercell atom label and cart. coord
        Rs = self.phonon.supercell.positions
        qs = self.get_commensurate_q_points_cartesian()
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
    def phonon(self) -> Phonopy:
        return self._phonon

    def set_at_qs(self, q_points: np.ndarray, fractional: bool = True):
        """prepare solution at given q_points via phonopy in ASE units

        Index notation:
            - i: (i, alpha) = primitive atom label and cart. coord
            - s: band index
            - q: q point index

        """
        if not fractional:  # convert to fractional coords for phonopy
            q_points = self.primitive.cell.reciprocal().scaled_positions(q_points)

        self.phonon.run_qpoints(
            q_points, with_eigenvectors=True, with_group_velocities=True,
        )
        data = self.phonon.get_qpoints_dict()

        self._set_data(data)
        self._q_points = q_points

    def _set_data(self, data: dict) -> None:
        """assign phonopy data to attributes"""
        w_sq = data["frequencies"].swapaxes(0, 1) / self.phonon._factor
        v_sq = data["group_velocities"].swapaxes(0, 1) / self.phonon._factor ** 2
        e_isq = np.moveaxis(data["eigenvectors"], 0, -1)

        if "weights" in data:
            weights_q = data["weights"]
        else:
            weights_q = np.ones(w_sq.shape[1])

        # make 0s 0, take 4th smallest absolute eigenvalue as reference
        thresh = sorted(abs(w_sq.flatten()))[3] * 0.1
        w_sq[abs(w_sq) < thresh] = 0
        w2_sq = np.sign(w_sq) * w_sq ** 2  # eigenvalues

        # build dynamical matrices
        ket_isq = e_isq[:, None, :, :]
        bra_jsq = e_isq[None, :, :, :].conj()
        D_ijq = (ket_isq * w2_sq[None, None, :] * bra_jsq).sum(axis=(2))

        # assign attributes
        self._w_sq, self._w2_sq, self._v_sq = w_sq, w2_sq, v_sq
        self._e_isq, self._D_ijq = e_isq, D_ijq
        self._weights_q = weights_q

        if (w_sq < 0).any():
            n_freq = len(w_sq[w_sq < 0])
            _talk(f"** {n_freq} negative frequencies! Check:")

    def set_at_commensurate_qs(self):
        """set mesh etc. at commensurate q-points"""
        self.set_at_qs(self.get_commensurate_q_points())

    @property
    def q_points(self):
        """return current set of q_points in fractional coords"""
        return self._q_points.copy()

    @property
    def q_points_cartesian(self):
        """return current set of q_points in Cartesian coords"""
        qs = self.q_points
        return self.primitive.cell.reciprocal().cartesian_positions(qs)

    def get_commensurate_q_points_cartesian(self):
        """return commensurate q-points in Cartesian coords [Nq, 3]"""
        qs = self.get_commensurate_q_points()
        return self.primitive.cell.reciprocal().cartesian_positions(qs)

    def get_commensurate_q_points(self):
        """return commensurate q-points in fractional coords [Nq, 3]"""
        return self._commensurate_q_points.copy()

    @property
    def qij(self):
        """return dynamical matrices at commensurate q-points [Nq, 3 * Np, 3 * Np]"""
        return np.moveaxis(self._D_ijq.copy(), -1, 0)

    @property
    def ijq(self):
        """return dynamical matrices at commensurate q-points [3 * Np, 3 * Np, Nq]"""
        return self._D_ijq.copy()

    @property
    def v_sq(self):
        """fractional group velocities in [3 * Np, Nq, 3]"""
        return self._v_sq

    @property
    def v_sq_cartesian(self):
        """Cartesian group velocities in [3 * Np, Nq, 3]"""
        return self.primitive.cell.cartesian_positions(self.v_sq)

    @property
    def w_sq(self):
        """eigenfrequencies in ASE units [3 * Np, Nq]"""
        return self._w_sq.copy()

    @property
    def w2_sq(self):
        """eigenvalues [3 * Np, Nq]"""
        return self._w2_sq.copy()

    @property
    def w_inv_sq(self):
        """inverse eigenfrequencies [3 * Np, Nq] respecting sign and zeros"""
        w = self.w_sq
        thresh = sorted(abs(w.flatten()))[3] * 0.9
        w[abs(w) < thresh] = 1e20
        w_inv_sq = 1 / w
        w_inv_sq[abs(w_inv_sq) < 1e-9] = 0
        return w_inv_sq

    @property
    def e_isq(self):
        """|i, sq> [3 * Np, 3 * Np, Nq]"""
        return self._e_isq.copy()

    @property
    def e_Isq(self):
        """eigenvectors w.r.t. comm. q-points |I, sq> [3 * Na, 3 * Np, Nq]"""
        return self._e_Isq.copy()

    @property
    def e_sqI(self):
        """project on eigenmodes <sq, I| [3 * Np, Nq, 3 * Na]"""
        return np.moveaxis(self.e_Isq, 0, -1)

    @property
    def weights_q(self):
        """return irreducible q-points weights [Nq_irred]"""
        return self._weights_q.copy()

    def set_mesh(self, mesh: np.ndarray, is_mesh_symmetry: bool = True):
        """set new q-points mesh via phonopy"""
        self.phonon.run_mesh(
            mesh,
            is_mesh_symmetry=is_mesh_symmetry,
            with_eigenvectors=True,
            with_group_velocities=True,
        )
        data = self.phonon.get_mesh_dict()
        self._q_points = data["qpoints"]
        self._set_data(data)

    def get_q_grid(self) -> tuple:
        """get q_points grid with symmetry information"""

        return get_ir_grid(self.q_points, primitive=self.primitive)

    def copy(self):
        """return copy"""
        return DynamicalMatrix(
            force_constants=self._fc_phonopy,
            primitive=self.primitive,
            supercell=self.supercell,
        )
