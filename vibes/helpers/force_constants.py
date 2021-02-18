""" Tools for dealing with force constants """

import numpy as np
from ase import Atoms
from ase.geometry import get_distances
from phonopy import Phonopy

# from vibes.ase.calculators.fc import get_smallest_vectors
from vibes.helpers import Timer, lazy_property, progressbar, talk
from vibes.helpers.lattice_points import get_lattice_points, map_I_to_iL
from vibes.helpers.numerics import clean_matrix
from vibes.helpers.supercell import map2prim
from vibes.helpers.supercell.supercell import supercell as fort
from vibes.io import get_identifier


_prefix = "force_constants"


def _talk(msg, verbose=True):
    return talk(msg, prefix=_prefix, verbose=verbose)


class ForceConstants:
    """Represents (phonopy) force constants and bundles related functionality"""

    def __init__(self, force_constants: np.ndarray, primitive: Atoms, supercell: Atoms):
        """instantiate ForceConstants with fc. in phonopy shape + ref. structures

        Args:
            force_constants: the force constant matrix in phonopy shape ([Np, Ns, 3, 3])
            primitive: the reference primitive structure
            supercell: the reference supercell

        """
        assert isinstance(supercell, Atoms)
        self.supercell = supercell.copy()
        # if primitive is None:  # create prim. cell if not given
        #     primitive = standardize_cell(supercell, to_primitive=True)
        assert isinstance(primitive, Atoms)
        self.primitive = primitive.copy()
        self._fc_phonopy = force_constants.copy()

        # map from supercell to primitive cell
        self._sc2pc = map2prim(primitive=self.primitive, supercell=self.supercell)

        # get lattice points in supercell
        lps, elps = get_lattice_points(primitive.cell, supercell.cell)
        self._lattice_points, self._lattice_points_extended = lps, elps

        # get mapping to lattice points
        I2iL, iL2I = map_I_to_iL(primitive, supercell, lattice_points=lps)
        self._I2iL, self._iL2I = I2iL, iL2I

        # get smallest vectors and their weights
        # self.smallest_vectors = get_smallest_vectors(self.supercell, self.primitive)

        # attach phonopy object
        smatrix = np.rint(supercell.cell @ primitive.cell.reciprocal().T).T
        phonon = Phonopy(unitcell=primitive, supercell_matrix=smatrix)
        phonon.force_constants = self._fc_phonopy
        self._phonon = phonon

    def __repr__(self):
        dct = get_identifier(self.primitive)
        sg, name = dct["space_group"], dct["material"]
        rep = f"Force Constants for {name} (sg {sg}) of shape {self._fc_phonopy.shape}"
        return repr(rep)

    @property
    def array(self) -> np.ndarray:
        """view on the raw phonopy-like force constants array"""
        return self._fc_phonopy.copy()

    @property
    def phonon(self) -> Phonopy:
        return self._phonon

    def get_remapped_to(self, atoms: Atoms, symmetrize: bool = True) -> np.ndarray:
        return remap_force_constants(
            force_constants=self._fc_phonopy,
            primitive=self.primitive,
            supercell=self.supercell,
            new_supercell=atoms,
            symmetrize=symmetrize,
        )

    @lazy_property
    def remapped_to_supercell(self) -> np.ndarray:
        return self.get_remapped_to(self.supercell)

    @property
    def remapped_to_supercell_3N_by_3N(self) -> np.ndarray:
        N, *_ = self.remapped_to_supercell.shape
        return self.remapped_to_supercell.swapaxes(1, 2).reshape((3 * N, 3 * N))

    @property
    def remapped(self) -> np.ndarray:
        """shorthand for remapped_to_supercell_3N_by_3N"""
        return self.remapped_to_supercell_3N_by_3N

    @property
    def reshaped_to_lattice_points(self) -> np.ndarray:
        return reshape_force_constants(
            primitive=self.primitive,
            supercell=self.supercell,
            force_constants=self.remapped_to_supercell,
            lattice_points=self.lattice_points,
        )

    @property
    def lattice_points(self) -> np.ndarray:
        return self._lattice_points

    @property
    def lattice_points_extended(self) -> list:
        return self._lattice_points_extended

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
        p2s_map = np.zeros(len(primitive), dtype=int)

        primitive.cell = new_supercell.cell

        new_supercell.wrap(eps=tol)
        primitive.wrap(eps=tol)

        for aa, a1 in enumerate(primitive):
            diff = new_supercell.positions - a1.position
            p2s_map[aa] = np.where(np.linalg.norm(diff, axis=1) < tol)[0][0]

        primitive.cell = primitive_cell
        primitive.wrap(eps=tol)

        return reduce_force_constants(fc_out, p2s_map)

    return fc_out


def reduce_force_constants(fc_full: np.ndarray, map2prim: np.ndarray):
    """reduce force constants from [N_sc, N_sc, 3, 3] to [N_prim, N_sc, 3, 3]

    Args:
        fc_full: The non-reduced force constant matrix
        map2prim: map from supercell to unitcell index

    Returns:
        The reduced force constants

    """
    uc_index = np.unique(map2prim)
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
        """see ForceConstants.__init__, but with mass weighting"""
        self.tol = tol
        # init ForceConstants
        super().__init__(
            force_constants=force_constants, primitive=primitive, supercell=supercell,
        )
        # now scale the forceconstants array
        w = primitive.get_masses()[:, None] * supercell.get_masses()[None, :]
        self._fc_phonopy /= w[:, :, None, None] ** 0.5

        # prepare fourier transforms q-space
        # prepare matrices
        n_lps = len(self.lattice_points)
        n_prim = len(self.primitive)
        w_Lm = np.zeros((n_lps, 8))  # weight of lattice point in extended supercell
        R_ijLm = np.zeros((n_prim, n_prim, n_lps, 8, 3))  # pair vectors incl. MIC imag.

        R_ij = self.primitive.positions[:, None] - self.primitive.positions[None, :]

        # get pair vectors which are within extended supercell
        cart2frac = self.supercell.cell.scaled_positions  # for fractional coords
        frac2cart = self.supercell.cell.cartesian_positions
        for LL, shell in enumerate(self.lattice_points_extended):
            for mm, lp in enumerate(shell):
                w_Lm[LL, mm] = 1 / len(shell)
                R_ijLm_frac = cart2frac(-R_ij.reshape(-1, 3) + lp).reshape(R_ij.shape)

                # map to to extended supercell [-0.5, 0.5] in frac. coords
                R1 = R_ijLm_frac
                d = 0.5 - tol * np.sign(R1)
                R = (R1 + d) % 1 - d
                R_ijLm[:, :, LL, mm] = frac2cart(R)

        # sanity check:
        shape = R_ijLm.shape
        R_ijLm_frac = cart2frac(R_ijLm.reshape(-1, 3)).reshape(shape)
        rmin, rmax = R_ijLm_frac.min(), R_ijLm_frac.max()
        assert rmin > -0.5 - tol, rmin
        assert rmax < +0.5 + tol, rmax

        # attach
        self._w_Lm = w_Lm  # [L, m]
        self._dm_ijL = self.reshaped_to_lattice_points[:, 0]  # [i, j, L, a, b]
        self._R_ijLm = R_ijLm  # [i, j, L, m, a]

    def at_qs(self, q_points: np.ndarray) -> np.ndarray:
        """compute dynamical matrices for a list of q-points"""
        _ = None  # shorthand for vector notation
        n_prim = len(self.primitive)
        pcell_inv = self.primitive.cell.reciprocal().T
        q_points_frac = (pcell_inv[_, :] * q_points[:, _, :]).sum(axis=-1)  # -> [q, a]

        # compute phases accounting for MIC images according to their weight
        e_qijLm = (q_points_frac[:, _, _, _, _, :] * self._R_ijLm[_, :]).sum(axis=-1)
        # -> [q, i, j, L, m]
        p_qijLm = np.exp(-2j * np.pi * e_qijLm)
        p_qijL = (self._w_Lm[_, _, _, :, :] * p_qijLm).sum(axis=-1)  # -> [q, i, j, L]

        Dq_L = p_qijL[:, :, :, :, _, _] * self._dm_ijL[_, :]  # [q, i, j, L a, b]
        Dq = Dq_L.sum(axis=3)  # -> [q, i, j, a, b]

        # reshape to [Nq, 3N, 3N]
        Dq = Dq.swapaxes(2, 3).reshape(len(q_points), 3 * n_prim, 3 * n_prim)

        # check if hermitian
        assert np.linalg.norm(Dq - Dq.swapaxes(1, 2).conj()) < self.tol

        return Dq

    def at_q(self, q_point: np.ndarray) -> np.ndarray:
        """return dynamical matrix at given q-point"""
        return self.at_qs(np.array([q_point]))[0]