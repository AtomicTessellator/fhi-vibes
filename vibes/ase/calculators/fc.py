import collections

import numpy as np
from ase import Atoms
from ase import units as u
from ase.calculators.calculator import Calculator
from ase.constraints import full_3x3_to_voigt_6_stress
from ase.geometry import find_mic
from phonopy.structure.cells import get_smallest_vectors as _get_smallest_vectors


def get_smallest_vectors(
    supercell: Atoms, primitive: Atoms, symprec: float = 1e-5
) -> collections.namedtuple:
    """get all pair vectors between atoms in supercell and primitive cell

    the routine respects m.i.c. and possible multiplicities when atoms are located
    at the boundary of the supercell.

    Returns:
        svecs, svecs_frac, multi, weights: smallest vectors in Cartesian and fractional
        coordinates (w.r.t. supercell), multiplicities, and weights (inverse mult.)

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
    svecs_frac, multi = _get_smallest_vectors(
        supercell_bases=supercell_bases,
        supercell_pos=supercell_pos,
        primitive_pos=primitive_pos,
        symprec=symprec,
    )

    # make transformation matrix integer and sanity check
    trans_mat_float = np.dot(supercell_bases, np.linalg.inv(primitive_bases))
    trans_mat = np.rint(trans_mat_float).astype(int)
    assert (np.abs(trans_mat_float - trans_mat) < 1e-8).all()
    svecs_frac = np.array(np.dot(svecs_frac, trans_mat), dtype="double", order="C")

    # compute Cartesian vectors
    svecs = supercell.cell.cartesian_positions(svecs_frac.reshape(-1, 3)).reshape(
        svecs_frac.shape
    )

    # compute weights based on multiplicity
    weights = np.zeros(svecs.shape[:3])
    for (ii, jj) in np.ndindex(svecs.shape[:2]):
        m = multi[ii, jj]
        weights[ii, jj, :m] = 1 / m

    return collections.namedtuple(
        "smallest_vectors", ("svecs", "svecs_frac", "multi", "weights")
    )(svecs, svecs_frac, multi, weights)


class FCCalculator(Calculator):
    """Harmonic Force Constants Calculator"""

    implemented_properties = ["energy", "energies", "forces"]
    implemented_properties += ["stress", "stresses", "heat_flux"]  # bulk properties

    def __init__(
        self,
        atoms_reference: Atoms = None,
        force_constants: np.ndarray = None,
        force_constants_qha: np.ndarray = None,
    ):
        r"""
        Args:
            atoms_reference: reference structure
            force_constants [3*N, 3*N]: force contants w.r.t. the reference structure
            force_constants_qha [3*N, 3*N]: quasiharmonic force constants d fc / d Vol

        Background:
            In 3D periodic systems, compute virial stresses according to

                s_i = -1/2V \sum_j (r_i - r_j) \otimes f_ij

            correctly accounting for periodic images and their multiplicity.
            Otherwise compute stresss according to

                s_i = -1/V u_i \otimes f_i

            which is sufficient in non-periodic systems.


                u_i: displacement of atom i
                f_i: force on atom i
                f_ij: force on atom i stemming from atom j

        """
        Calculator.__init__(self)

        self.atoms_reference = atoms_reference
        self.force_constants = force_constants
        self.pbc = self.atoms_reference.cell.rank == 3
        # reshape force constants from [Ia, Jb] -> [I, J, a, b]
        fc_shape = (len(atoms_reference), 3, len(atoms_reference), 3)
        self.fc_IJab = force_constants.reshape(*fc_shape).swapaxes(1, 2)

        if force_constants_qha is not None:
            self.qha = True
            self.fc_qha_IJab = force_constants_qha.reshape(*fc_shape).swapaxes(1, 2)
        else:
            self.qha = False

        if self.pbc:
            # save pair vectors and multiplicities
            atoms = atoms_reference
            self.pairs, *_, self.weights = get_smallest_vectors(atoms, atoms)

            # pair vectors w.r.t. reference pos. including multiplicites
            r0_ijm = self.pairs  # [I, J, m, a]
            # average over equivalent periodic images m
            self.r0_ij = (self.weights[:, :, :, None] * r0_ijm).sum(
                axis=2
            )  # -> [I, J, a]

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
        if self.pbc:
            volume = atoms.get_volume()
        else:
            volume = 1.0

        # prepare calculation using following notation for indices:
        # I, J: 1 ... N (atom index)
        # a, b: 1 ... 3 (vector component)
        # m: 1 ... M (m.i.c. multiplicity, at most 27)
        dr = displacements  # [I, a]
        # reshaped force constants [I, J, a, b]
        fc = self.fc_IJab
        # PU = sum_b forceconstants P [I, J, a, b] * displacements U [J, b]
        PU = (fc * dr[None, :, None, :]).sum(axis=-1)  # -> [I, J, a]

        # i) account for m.i.c. multiplicity
        if self.pbc:
            # pair displacements as matrix including multiplicites
            d_ijm = (dr[:, None, :] - dr[None, :, :])[:, :, :, None]  # [I, J, a, m]
            # average over equivalent periodic images m (mask comes w/ factor 1/M)
            d_ij = (self.weights[:, :, None, :] * d_ijm).sum(axis=-1)  # -> [I, J, a]

            # stress matrix s_ij = r_ij [I, J, a] * PU [I, J, b] / volume
            # stress term from reference positions for heat flux:
            s0_ij = self.r0_ij[:, :, :, None] * PU[:, :, None, :]  # -> [I, J, a, b]
            s0_ij /= volume
            # stress term due to displacements for stress/pressure
            sd_ij = d_ij[:, :, :, None] * PU[:, :, None, :]  # -> [I, J, a, b]
            sd_ij /= volume

        # ii) w/o p.b.c or in very simple force fields life is easy:
        else:
            # reshape force constants from 3Nx3N to NxNx3x3
            # s_ij [I, J, a, b] = u_i [I, :, a, :] * PU [I, J, :, b]
            sd_ij = dr[:, None, :, None] * PU[:, :, None, :]  # -> [I, J, a, b]
            sd_ij /= volume
            s0_ij = np.zeros_like(sd_ij)

        # from here on the two cases are similar
        # QHA contribution (similar to ha. case, but prefactor 3/2 instead of 1/V)
        if self.qha:
            fc_qha = self.fc_qha_IJab  # [I, J, a, b]
            PU_qha = (fc_qha * dr[None, :, None, :]).sum(axis=-1)  # [I, J, a]

            s_ij_qha = dr[:, None, :, None] * PU_qha[:, :, None, :]  # -> [I, J, a, b]
            s_ij_qha *= 3 / 2

            # add to displacemnt-stress sd:
            sd_ij += s_ij_qha

        # assign properties
        # potential energy [I] = 1/2 * sum_a - dr [I, a] * f [I, a]
        energies = -(displacements * forces).sum(axis=-1) / 2  # -> [I]
        energy = energies.sum()  # -> [1]

        self.results["forces"] = forces
        self.results["energy"] = energy
        self.results["energies"] = energies

        # no lattice, no stress
        if self.pbc:
            #  atomic stresses s_i [I, a, b] = sum_j s_ij [I, J, a, b]
            s0_i = s0_ij.sum(axis=1)  # -> [I, a, b]
            sd_i = sd_ij.sum(axis=1)  # -> [I, a, b]
            s_i = s0_i + sd_i

            stress = full_3x3_to_voigt_6_stress(s_i.sum(axis=0))
            self.results["stress"] = stress
            self.results["stresses"] = s_i
            self.results["stresses_0"] = s0_i

            if self.qha:
                self.results["stresses_qha"] = s_ij_qha.sum(axis=1)

            # # compute heat flux in eV / AA**2 / fs
            # J = 1/2 * sum_i s_i \cdot p_i
            pref = 0.5
            v = atoms.get_velocities() * u.fs
            j = pref * (s_i * v[:, None, :]).sum(axis=(0, -1))
            j0 = pref * (s0_i * v[:, None, :]).sum(axis=(0, -1))
            self.results["heat_flux"] = j
            self.results["heat_flux_0"] = j0
        return self.results