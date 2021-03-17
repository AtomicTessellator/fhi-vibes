""" Tools for dealing with force constants """

import numpy as np
from ase import Atoms
from ase.geometry import get_distances
from phonopy import Phonopy

from vibes.brillouin import get_ir_grid
from vibes.helpers import Timer, lazy_property, progressbar, talk
from vibes.helpers.lattice_points import (
    get_commensurate_q_points,
    get_lattice_points,
    map_I_to_iL,
)
from vibes.helpers.numerics import clean_matrix
from vibes.helpers.supercell import map2prim
from vibes.helpers.supercell.supercell import supercell as fort
from vibes.io import get_identifier


la = np.linalg


_prefix = "force_constants"


def _talk(msg, verbose=True):
    return talk(msg, prefix=_prefix, verbose=verbose)


class ForceConstants:
    """Represents (phonopy) force constants and bundles related functionality"""

    def __init__(self, force_constants: np.ndarray, primitive: Atoms, supercell: Atoms):
        """instantiate ForceConstants with fc. in phonopy shape + ref. structures

        Args:
            force_constants: the force constant matrix in phonopy shape ([Np, Ns, 3, 3])
                             or in remapped shape [3 * Na, 3 * Na]
            primitive: the reference primitive structure
            supercell: the reference supercell

        """
        assert isinstance(supercell, Atoms)
        self.supercell = supercell.copy()
        # if primitive is None:  # create prim. cell if not given
        #     primitive = standardize_cell(supercell, to_primitive=True)
        assert isinstance(primitive, Atoms)
        self.primitive = primitive.copy()

        # map from supercell to primitive cell
        self._sc2pc = map2prim(primitive=self.primitive, supercell=self.supercell)

        Np, Na = len(primitive), len(supercell)
        if force_constants.shape == (3 * Na, 3 * Na):
            fc_full = force_constants.reshape([Na, 3, Na, 3]).swapaxes(1, 2)
            force_constants = reduce_force_constants(fc_full, map2prim=self._sc2pc)
        elif force_constants.shape == (Na, Na, 3, 3):
            fc_full = force_constants
            force_constants = reduce_force_constants(fc_full, map2prim=self._sc2pc)
        else:
            assert np.shape(force_constants) == (Np, Na, 3, 3)

        self._fc_phonopy = force_constants.copy()
        self._mass_weights = np.ones([len(primitive), len(supercell)])

        # get lattice points in supercell
        lps, elps = get_lattice_points(primitive.cell, supercell.cell)
        self._lattice_points, self._lattice_points_extended = lps, elps

        # get mapping to lattice points
        I2iL, iL2I = map_I_to_iL(primitive, supercell, lattice_points=lps)
        self._I2iL, self._iL2I = I2iL, iL2I

    def __repr__(self):
        dct = get_identifier(self.primitive)
        sg, name = dct["space_group"], dct["material"]
        rep = f"Force Constants for {name} (sg {sg}) of shape {self.array.shape}"
        return repr(rep)

    # data
    @property
    def array(self) -> np.ndarray:
        """view on the raw phonopy-like force constants array INCLUDING mass-weighting"""
        return self._fc_phonopy.copy() * self._mass_weights[:, :, None, None]

    # remapping to supercell
    @lazy_property
    def remapped_to_supercell(self) -> np.ndarray:
        """"return in [I, J, a, b] including symmetrization"""
        return remap_force_constants(
            force_constants=self.array,
            primitive=self.primitive,
            supercell=self.supercell,
            new_supercell=self.supercell,
            symmetrize=True,
        )

    @property
    def remapped_to_supercell_3N_by_3N(self) -> np.ndarray:
        """"return in [I, a, J, b]"""
        N, *_ = self.remapped_to_supercell.shape
        return self.remapped_to_supercell.swapaxes(1, 2).reshape((3 * N, 3 * N))

    @property
    def remapped(self) -> np.ndarray:
        """"shorthand for `remapped_to_supercell_3N_by_3N, return in [I, a, J, b]"""
        return self.remapped_to_supercell_3N_by_3N

    # lattice points things
    @property
    def I2iL_map(self):
        return self._I2iL.copy()

    @property
    def lattice_points(self) -> np.ndarray:
        """return list of lattice points contained in supercell"""
        return self._lattice_points

    @lazy_property
    def reshaped_to_lattice_points(self) -> np.ndarray:
        """"return in [i, L, j, K, a, b"""
        return reshape_force_constants(
            primitive=self.primitive,
            supercell=self.supercell,
            force_constants=self.remapped_to_supercell,
            lattice_points=self.lattice_points,
        )

    def get_dynamical_matrix(self) -> object:
        """return mass-weighted forceconstants as DynamicalMatrix object"""
        return DynamicalMatrix(
            force_constants=self._fc_phonopy,
            primitive=self.primitive,
            supercell=self.supercell,
        )


def remap_force_constants(
    force_constants: np.ndarray,
    primitive: Atoms,
    supercell: Atoms,
    new_supercell: Atoms = None,
    reduce_fc: bool = False,
    two_dim: bool = False,
    symmetrize: bool = True,
    tol: float = 1e-5,
    eps: float = 1e-13,
    fortran: bool = True,
) -> np.ndarray:
    """remap force constants [N_prim, N_sc, 3, 3] to [N_sc, N_sc, 3, 3]

    Args:
        force_constants: force constants in [N_prim, N_sc, 3, 3] shape
        primitive: primitive cell for reference
        supercell: supercell for reference
        new_supercell: supercell to map to
        reduce_fc: return in [N_prim, N_sc, 3, 3]  shape
        two_dim: return in [3*N_sc, 3*N_sc] shape
        symmetrize: make force constants symmetric
        tol: tolerance to discern pairs
        eps: finite zero

    Returns:
        The remapped force constants

    """
    timer = Timer("remap force constants", prefix=_prefix)

    if new_supercell is None:
        new_supercell = supercell.copy()

    # find positions of primitive cell within the (reference) supercell
    primitive_cell = primitive.cell.copy()
    primitive.cell = supercell.cell

    primitive.wrap(eps=tol)
    supercell.wrap(eps=tol)

    n_sc_new = len(new_supercell)

    # make a list of all pairs for each atom in primitive cell
    sc_r = np.zeros((force_constants.shape[0], force_constants.shape[1], 3))
    for aa, a1 in enumerate(primitive):
        diff = supercell.positions - a1.position
        p2s = np.where(np.linalg.norm(diff, axis=1) < tol)[0][0]
        # sc_r[aa] = supercell.get_distances(p2s, range(n_sc), mic=True, vector=True)
        # replace with post 3.18 ase routine:
        spos = supercell.positions
        sc_r[aa], _ = get_distances([spos[p2s]], spos, cell=supercell.cell, pbc=True)

    # find mapping from supercell to origin atom in primitive cell
    primitive.cell = primitive_cell
    sc2pc = map2prim(primitive, new_supercell)

    if fortran:
        fc_out = fort.remap_force_constants(
            positions=new_supercell.positions,
            pairs=sc_r,
            fc_in=force_constants,
            map2prim=sc2pc,
            inv_lattice=new_supercell.cell.reciprocal(),
            tol=tol,
            eps=eps,
        )
    else:

        _talk(".. use python")
        ref_struct_pos = new_supercell.get_scaled_positions(wrap=True)
        sc_temp = new_supercell.get_cell(complete=True)

        fc_out = np.zeros((n_sc_new, n_sc_new, 3, 3))
        for a1, r0 in enumerate(progressbar(new_supercell.positions)):
            uc_index = sc2pc[a1]
            for sc_a2, sc_r2 in enumerate(sc_r[uc_index]):
                r_pair = r0 + sc_r2
                r_pair = np.linalg.solve(sc_temp.T, r_pair.T).T % 1.0
                for a2 in range(n_sc_new):
                    r_diff = np.abs(r_pair - ref_struct_pos[a2])
                    # Integer value is the equivalent of 0.0
                    r_diff -= np.floor(r_diff + eps)
                    if np.linalg.norm(r_diff) < tol:
                        fc_out[a1, a2, :, :] += force_constants[uc_index, sc_a2, :, :]

    timer()

    # check symmetries:
    fc_2d = fc_out.swapaxes(1, 2).reshape((3 * n_sc_new, 3 * n_sc_new))  # -> 3Nx3N
    violation = np.linalg.norm(fc_2d - fc_2d.T)
    if violation > 1e-5:
        _talk(f"** Force constants are not symmetric by {violation:.2e}.")
    if symmetrize:
        _talk("-> Symmetrize force constants.")
        fc_2d = 0.5 * (fc_2d + fc_2d.T)

    # sum rule 1
    violation = abs(fc_2d.sum(axis=0)).mean()
    if violation > 1e-9:
        _talk(f"** Sum rule violated by {violation:.2e} (axis 1).")

    if two_dim:
        return fc_2d

    fc_out = fc_2d.reshape((n_sc_new, 3, n_sc_new, 3)).swapaxes(1, 2)  # -> NxNx3x3

    if reduce_fc:
        s2p_map = map2prim(primitive, new_supercell)

        return reduce_force_constants(fc_out, s2p_map)

    return fc_out


def reduce_force_constants(fc_full: np.ndarray, map2prim: np.ndarray):
    """reduce force constants from [N_sc, N_sc, 3, 3] to [N_prim, N_sc, 3, 3]

    Args:
        fc_full: The non-reduced force constant matrix
        map2prim: map from supercell to unitcell index

    Returns:
        The reduced force constants

    """
    _, uc_index = np.unique(map2prim, return_index=True)
    fc_out = np.zeros((len(uc_index), fc_full.shape[1], 3, 3))
    for ii, uc_ind in enumerate(uc_index):
        fc_out[ii, :, :, :] = fc_full[uc_ind, :, :, :]

    return fc_out


def reshape_force_constants(
    primitive: Atoms,
    supercell: Atoms,
    force_constants: np.ndarray,
    scale_mass: bool = False,
    lattice_points: np.ndarray = None,
    symmetrize: bool = True,
) -> np.ndarray:
    """ reshape from (N_prim x N_super x 3 x 3) into 3x3 blocks labelled by (i,L)

    Args:
        primitive: The primitive cell structure
        supercell: The super cell structure
        force_constants: The input force constant matrix
        scale_mass: If True scale phi by the product of the masses
        lattice_points: the lattice points to include

    Returns:
        The remapped force constants

    """
    if len(force_constants.shape) > 2:
        if force_constants.shape[0] != force_constants.shape[1]:
            force_constants = remap_force_constants(
                force_constants,
                primitive,
                supercell,
                two_dim=True,
                symmetrize=symmetrize,
            )
        else:
            n_sc = force_constants.shape[0]
            force_constants = force_constants.swapaxes(1, 2).reshape(
                (3 * n_sc, 3 * n_sc)
            )

    if lattice_points is None:
        lattice_points, _ = get_lattice_points(primitive.cell, supercell.cell)

    indeces, _ = map_I_to_iL(primitive, supercell, lattice_points=lattice_points)

    n_i = len(primitive)
    n_L = len(lattice_points)

    masses = primitive.get_masses()

    new_force_constants = np.zeros([n_i, n_L, n_i, n_L, 3, 3])

    for n1 in range(len(supercell)):
        for n2 in range(len(supercell)):
            phi = force_constants[3 * n1 : 3 * n1 + 3, 3 * n2 : 3 * n2 + 3]

            i1, L1, i2, L2 = (*indeces[n1], *indeces[n2])

            if scale_mass:
                phi /= np.sqrt(masses[i1] * masses[i2])

            new_force_constants[i1, L1, i2, L2] = clean_matrix(phi)

    return new_force_constants


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