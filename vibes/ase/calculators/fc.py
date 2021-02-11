import numpy as np
from ase import Atoms
from ase import units as u
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
        force_constants_qha: np.ndarray = None,
        virials: bool = True,
    ):
        r"""
        Args:
            atoms_reference: reference structure
            force_constants [3*N, 3*N]: force contants w.r.t. the reference structure
            force_constants_qha [3*N, 3*N]: quasiharmonic force constants d fc / d Vol
            virials: compute stress via virial formula

        Background:
            When `virials` is chosen, compute virial stresses according to

                s_i = -1/2V \sum_j (r_i - r_j) \otimes f_ij

            correctly accounting for periodic images and their multiplicity.
            Otherwise compute stresss according to

                s_i = -1/V u_i \otimes f_i

            which discards the relative position of atoms. This might be besser to
            approximate pressure, but the stress tensor components will be off.

            where
                u_i: displacement of atom i
                f_i: force on atom i
                f_ij: force on atom i stemming from atom j

        """
        Calculator.__init__(self)

        self._virials = virials
        self.atoms_reference = atoms_reference
        self.force_constants = force_constants

        if force_constants_qha is not None:
            self.qha = True
            self.force_constants_qha = force_constants_qha
        else:
            self.qha = False
            self.force_constants_qha = np.zeros_like(force_constants)

        if self._virials:
            # save pair vectors and multiplicities
            atoms = atoms_reference
            self.pairs, self.mult = get_smallest_vectors(atoms, atoms)

            # create a mask to correctly account for multiple images
            self.mask = np.zeros(self.pairs.shape[:3])
            for (ii, jj) in np.ndindex(self.pairs.shape[:2]):
                m = self.mult[ii, jj]
                self.mask[ii, jj, :m] = 1 / m

    def calculate(self, atoms, *unused_args, **unused_kw):

        n_atoms = len(self.atoms_reference)

        assert n_atoms == len(atoms), f"Atom number changed ({n_atoms} -> {len(atoms)})"

        # find displacements respecting p.b.c. and m.i.c.
        cell = np.asarray(atoms.cell)
        positions = atoms.positions
        reference_positions = self.atoms_reference.positions
        displacements = find_mic(positions - reference_positions, cell)[0]

        # forces are always f_i = - phi_ij @ u_j
        fc = self.force_constants
        forces = -(fc @ displacements.flatten()).reshape(displacements.shape)

        # volume
        if atoms.cell.rank == 3:
            volume = atoms.get_volume()
        else:
            volume = 1.0

        # prepare calculation using following notation for indices:
        # I, J: 1 ... N (atom index)
        # a, b: 1 ... 3 (vector component)
        # m: 1 ... M (m.i.c. multiplicity, at most 27)
        dr = displacements  # [I, a]
        # reshape force constants from [Ia, Jb] -> [I, J, a, b]
        fc_shape = (n_atoms, 3, n_atoms, 3)
        fc = self.force_constants.reshape(*fc_shape).swapaxes(1, 2)
        # PU = sum_b forceconstants P [I, J, a, b] * displacements U [J, b]
        PU = (fc * dr[None, :, None, :]).sum(axis=-1)  # -> [I, J, a]

        # i) account for m.i.c. multiplicity
        if self._virials:
            # pair displacements as matrix
            d_ij = dr[:, None, :] - dr[None, :, :]  # [I, J, a]
            # pair vectors w.r.t. reference pos. including multiplicites
            r0_ijm = self.pairs @ self.atoms_reference.cell
            # pair vectors including displacements and multiplicites
            r_ijm = r0_ijm + d_ij[:, :, None, :]  # [I, J, m, a]
            # average over equivalent periodic images m (mask comes w/ factor 1/M)
            r_ij = (self.mask[:, :, :, None] * r_ijm).sum(axis=2)  # -> [I, J, a]
            # dev:
            # r0_ij = (self.mask[:, :, :, None] * r0_ijm).sum(axis=2)  # -> [I, J, a]
            # dr_ij = (self.mask[:, :, :, None] * d_ij[:, :, None, :]).sum(axis=2)

            # stress matrix s_ij = r_ij [I, J, a] * PU [I, J, b] / volume
            s_ij = r_ij[:, :, :, None] * PU[:, :, None, :]  # -> [I, J, a, b]
            s_ij /= volume
            # dev:
            # s0_ij = r0_ij[:, :, :, None] * PU[:, :, None, :]  # -> [I, J, a, b]
            # sd_ij = dr_ij[:, :, :, None] * PU[:, :, None, :]  # -> [I, J, a, b]
            # s0_ij /= volume
            # sd_ij /= volume

        # ii) w/o p.b.c or in very simple force fields life is easy:
        else:
            # reshape force constants from 3Nx3N to NxNx3x3
            # s_ij [I, J, a, b] = u_i [I, :, a, :] * PU [I, J, :, b]
            s_ij = dr[:, None, :, None] * PU[:, :, None, :]  # -> [I, J, a, b]
            s_ij /= volume

        # from here on the two cases are similar
        # QHA contribution (similar to ha. case, but prefactor 3/2 instead of 1/V)
        if self.qha:
            fc_qha = self.force_constants_qha
            fc_qha = fc_qha.reshape(*fc_shape).swapaxes(1, 2)  # [I, J, a, b]
            PU_qha = (fc_qha * dr[None, :, None, :]).sum(axis=-1)  # [I, J, a]

            s_ij_qha = dr[:, None, :, None] * PU_qha[:, :, None, :]  # -> [I, J, a, b]
            s_ij_qha *= 3 / 2

        # assign properties
        # potential energy [I] = 1/2 * sum_a - dr [I, a] * f [I, a]
        energies = -(displacements * forces).sum(axis=1) / 2  # -> [I]
        energy = energies.sum()  # -> [1]

        self.results["forces"] = forces
        self.results["energy"] = energy
        self.results["energies"] = energies
        self.results["free_energy"] = energy

        # no lattice, no stress
        if atoms.cell.rank == 3:
            #  atomic stresses s_i [I, a, b] = sum_j s_ij [I, J, a, b]
            if self.qha:
                s_ij += s_ij_qha
            s_i = s_ij.sum(axis=1)  # -> [I, a, b]
            # dev
            # s0_i = s0_ij.sum(axis=1)  # -> [I, a, b]
            # sd_i = sd_ij.sum(axis=1)

            stress = full_3x3_to_voigt_6_stress(s_i.sum(axis=0))
            self.results["stress"] = stress
            self.results["stresses"] = s_i
            # dev
            # self.results["stresses_0"] = s0_i
            # self.results["stresses_d"] = sd_i
            if self.qha:
                self.results["stresses_qha"] = s_ij_qha.sum(axis=1)

            # # compute heat flux in eV / AA**2 / fs
            # J = 1/2 * sum_i s_i \cdot p_i
            pref = 0.5
            v = atoms.get_velocities() / u.fs
            j = pref * (s_i * v[:, None, :]).sum(axis=-1).sum(axis=0)
            self.results["heat_flux"] = j
