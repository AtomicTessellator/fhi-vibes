import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.constraints import full_3x3_to_voigt_6_stress
from ase.geometry import find_mic
from phonopy.structure.cells import get_smallest_vectors as _get_smallest_vectors


def get_smallest_vectors(supercell: Atoms, primitive: Atoms, symprec: float = 1e-5):
    """get all pair vectors between atoms in supercell and primitive cell

    the routine respects m.i.c. and possible multiplicities when atoms are located
    at the boundary of the supercell.

    Returns:
        svecs, mult: smallest vectors in fractional coordinates (w.r.t. supercell) and
            multiplicities

    """
    # save bases
    supercell_bases = supercell.cell
    primitive_bases = primitive.cell

    # get fractional positions w.r.t to supercell basis
    supercell_pos = supercell.get_scaled_positions()
    primitive_w_supercell_basis = primitive.copy()
    primitive_w_supercell_basis.cell = supercell.cell.copy()
    primitive_pos = primitive_w_supercell_basis.get_scaled_positions()

    # get the smallest vectors and multiplicites
    svecs, multi = _get_smallest_vectors(
        supercell_bases=supercell_bases,
        supercell_pos=supercell_pos,
        primitive_pos=primitive_pos,
        symprec=symprec,
    )

    # make transformation matrix integer and sanity check
    trans_mat_float = np.dot(supercell_bases, np.linalg.inv(primitive_bases))
    trans_mat = np.rint(trans_mat_float).astype(int)
    assert (np.abs(trans_mat_float - trans_mat) < 1e-8).all()

    # return
    svecs = np.array(np.dot(svecs, trans_mat), dtype="double", order="C")
    return svecs, multi


class FCCalculator(Calculator):
    """Harmonic Force Constants Calculator"""

    implemented_properties = ["energy", "energies", "forces", "free_energy"]
    implemented_properties += ["stress", "stresses"]  # bulk properties

    def __init__(
        self,
        atoms_reference: Atoms = None,
        force_constants: np.ndarray = None,
        pbc: bool = True,
        fast: bool = True,
        **kwargs,
    ):
        r"""
        Args:
            atoms_reference: reference structure
            force_constants [3*N, 3*N]: force contants w.r.t. the reference structure
            pbc: compute stress correctly accounting for p.b.c.
            fast: use vectorized algorithm to speed up computation of stress

        Background:
            When `pbc` is chosen, compute virial stresses according to

                s_i = -1/2 \sum_j (u_i - u_j) \otimes f_ij

            correctly accounting for periodic images and their multiplicity.
            Otherwise compute stresss according to

                s_i = -u_i \otimes f_i

            which is also correct for LennardJones-type force fields.

            where
                u_i: displacement of atom i
                f_i: force on atom i
                f_ij: force on atom i stemming from atom j

        """
        Calculator.__init__(self, **kwargs)

        self._pbc = pbc
        self._fast = fast
        self.atoms_reference = atoms_reference
        self.force_constants = force_constants

        if self._pbc:
            # save pair vectors and multiplicities
            atoms = atoms_reference
            self.pairs, self.mult = get_smallest_vectors(atoms, atoms)

            # create a mask to correctly account for multiple images
            self.mask = np.zeros(self.pairs.shape[:3])
            for (ii, jj) in np.ndindex(self.pairs.shape[:2]):
                m = self.mult[ii, jj]
                self.mask[ii, jj, :m] = 1

    def calculate(self, atoms, *unused_args, **unused_kw):

        n_atoms = len(self.atoms_reference)

        assert n_atoms == len(atoms), f"Atom number changed ({n_atoms} -> {len(atoms)})"

        # find displacements respecting p.b.c. and m.i.c.
        cell = np.asarray(atoms.cell)
        positions = atoms.positions
        reference_positions = self.atoms_reference.positions
        displacements = find_mic(positions - reference_positions, cell)[0]

        # forces are always f_i = - phi_ij @ dr_j
        fc = self.force_constants
        forces = -(fc @ displacements.flatten()).reshape(displacements.shape)

        # stress:es
        # i) vectorized
        if self._pbc and self._fast:
            # reshape force constants from 3Nx3N to NxNx3x3
            fc = self.force_constants.reshape(n_atoms, 3, n_atoms, 3).swapaxes(1, 2)
            dr = displacements[np.newaxis, :, np.newaxis, :]
            # forceconstants * displacements [N, N, 3]
            Du = (fc * dr).sum(axis=-1)
            # pair displacemnts as matrix [N, N, 3]
            d_ij = displacements[:, np.newaxis] - displacements[np.newaxis, :]
            # pair vectors w.r.t. reference pos. including multiplicites [N, N, 27, 3]
            r0_ijm = self.pairs @ self.atoms_reference.cell
            # pair vectors including displacements and multiplicites [N, N, 27, 3]
            r_ijm = r0_ijm + d_ij[:, :, None, :]
            # sum over equivalent periodic images [N, N, 3]
            r_ij = (self.mask[:, :, :, None] * r_ijm).sum(axis=2)
            # s_ij = r_ij * D_ij @ u_j [N, N, 3, 3]
            s_ij = r_ij[:, :, :, None] * Du[:, :, None, :]

            # s_i [N, 3, 3]
            stresses = s_ij.sum(axis=1)

            # dev: perform sum over multiple images later, more memory consuming
            # s_ijm = r_ijm[:, :, :, :, None] * Du[:, :, None, None, :]
            # s_ij = (self.mask[:, :, :, None, None] * s_ijm).sum(axis=2)
        # ii) naive implementation w/ loops
        elif self._pbc and not self._fast:
            forces = np.zeros((n_atoms, n_atoms, 3))
            stresses = np.zeros((n_atoms, n_atoms, 3, 3))
            for (ii, jj) in np.ndindex(n_atoms, n_atoms):
                n_mult = self.mult[ii, jj]
                d_ij = displacements[ii] - displacements[jj]
                p_ij = fc[3 * ii : 3 * (ii + 1), 3 * jj : 3 * (jj + 1)]
                f_ij = -p_ij @ displacements[jj]

                # REM: fully antisymmetric force f_ij = -f_ji looks like this
                # p_ji = fc[3 * jj : 3 * (jj + 1), 3 * ii : 3 * (ii + 1)]
                # f_ij = -p_ij @ displacements[jj] + p_ji @ displacements[ii]

                forces[ii, jj] = f_ij

                # sum up including equivalent images
                for mm in range(n_mult):
                    r_ij_m = (self.pairs[ii, jj, mm] @ self.atoms_reference.cell) + d_ij
                    stresses[ii, jj] -= np.outer(r_ij_m, f_ij)  # / 2 when antisym. f_ij

                # sanity check
                if ii == jj:
                    assert np.allclose(stresses[ii, jj], 0)

            forces = forces.sum(axis=1)
            stresses = stresses.sum(axis=1)

        # iii) w/o p.b.c or in very simple force fields life is easy:
        else:
            stresses = np.zeros((n_atoms, 3, 3))
            for ii in range(n_atoms):
                stresses[ii] = -np.outer(displacements[ii], forces[ii])

        # assign properties
        # potential energy = - dr * f
        energies = -(displacements * forces).sum(axis=1)
        energy = energies.sum()

        self.results["energy"] = energy
        self.results["energies"] = energies
        self.results["free_energy"] = energy

        self.results["forces"] = forces

        # no lattice, no stress
        if atoms.cell.rank == 3:
            stresses = full_3x3_to_voigt_6_stress(stresses)
            self.results["stress"] = stresses.sum(axis=0) / atoms.get_volume()
            self.results["stresses"] = stresses / atoms.get_volume()
