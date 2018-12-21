""" helps to find lattice points in supercell, match positions to images in the 
unit cell etc. """

from itertools import product
import numpy as np
import scipy.linalg as la
from hilde.helpers.timer import Timer
from hilde.helpers.maths import clean_matrix
from . import supercell as sc


def get_dynamical_matrix(q, primitive, supercell, force_constants, eps=1e-12):
    """ build the dynamical matrix for one q_point """
    return get_dynamical_matrices([q], primitive, supercell, force_constants, eps)[0]


def _prefactor(q, l):
    return np.exp(-2j * np.pi * q @ l)


def get_dynamical_matrices(q_points, primitive, supercell, force_constants, eps=1e-12):
    """ build the dynamical matrix for each q_point """

    lattice_points, lattice_points_ext = get_lattice_points(primitive, supercell)

    force_constants_reshaped = reshape_force_constants(
        primitive, supercell, force_constants, lattice_points=lattice_points
    )

    masses = primitive.get_masses()

    n_prim = len(primitive)

    dyn_matrices = []

    for qq in q_points:
        dyn_matrix = np.zeros([n_prim, n_prim, 3, 3], dtype=complex)

        for LL, lp in enumerate(lattice_points):
            prefactor = _prefactor(qq, lp)
            if np.linalg.norm(prefactor.imag) < eps:
                prefactor = prefactor.real

            for ii, ri in enumerate(primitive.positions):
                for jj, rj in enumerate(primitive.positions):
                    dyn_matrix[ii, jj, :, :] += (
                        0.5 * prefactor * force_constants_reshaped[ii, 0, jj, LL]
                    )

                    dyn_matrix[ii, jj, :, :] += (
                        0.5 * prefactor * force_constants_reshaped[jj, LL, ii, 0]
                    )

        for ii in range(n_prim):
            for jj in range(n_prim):
                dyn_matrix[ii, jj, :, :] /= np.sqrt(masses[ii] * masses[jj])

        dyn_matrices.append(dyn_matrix.swapaxes(1, 2).reshape(3 * n_prim, 3 * n_prim))

    return dyn_matrices


def reshape_force_constants(
    primitive, supercell, force_constants, scale_mass=False, lattice_points=None
):
    """ reshape from (3N x 3N) into 3x3 blocks labelled by (i,L) """

    if lattice_points is None:
        lattice_points, _ = get_lattice_points(primitive, supercell)

    indeces = assign_primitive_positions(
        primitive, supercell, lattice_points=lattice_points
    )

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


def get_commensurate_q_points(atoms, supercell, tolerance=1e-5, **kwargs):
    """ For a commensurate q_points we have

        exp( 2\pi q . L_k ) = 1 for any k and L_k being the supercell lattice vectors

        in other workds, q is a linear combination of G_k, where G_k are the inverse
        lattice vectors of the supercell lattice. Only thos are counted which fit into
        the inverse lattice of the primitive cell.
        This means we have to call lattice_points.get_lattice_points with the inverse
        lattices.

    """

    lattice = atoms.cell.T
    superlattice = supercell.cell.T

    inv_lattice = la.inv(lattice)
    inv_superlattice = la.inv(superlattice)

    inv_lattice_points, _ = _get_lattice_points(inv_superlattice, inv_lattice, **kwargs)

    return inv_lattice_points


def assign_primitive_positions(
    in_atoms, in_supercell, lattice_points=None, extended=False, tolerance=1e-5
):
    """ write positions in the supercell as (i, L), where i is the respective index
        in the primitive cell and L is the lattice point """

    timer = Timer()

    atoms = in_atoms.copy()
    supercell = in_supercell.copy()
    atoms.wrap()
    supercell.wrap()

    if lattice_points is None:
        if extended:
            _, lattice_points = get_lattice_points(atoms, supercell)
        else:
            lattice_points, _ = get_lattice_points(atoms, supercell)

    # create all positions R = r_i + L
    all_positions = []
    tuples = []
    for ii, pos in enumerate(atoms.positions):
        for LL, lp in enumerate(lattice_points):
            all_positions.append(pos + lp)
            tuples.append((ii, LL))

    indices = []
    matches = []
    for jj, spos in enumerate(supercell.positions):
        for ii, pos in enumerate(atoms.positions):
            for LL, lp in enumerate(lattice_points):
                if la.norm(spos - pos - lp) < tolerance:
                    indices.append((ii, LL))
                    matches.append(jj)
                    break

    # catch possibly unwrapped atoms
    for jj, spos in enumerate(supercell.positions):
        if jj in matches:
            continue
        for LL, lp in enumerate(lattice_points):
            for ii, pos in enumerate(atoms.positions):
                fpos = fractional(spos - pos - lp, supercell.cell)
                tol = tolerance
                if la.norm((fpos + tol) % 1 - tol) < tolerance:
                    indices.append((ii, LL))
                    matches.append(jj)
                    break

    # sanity checks:
    if len(np.unique(matches)) != len(supercell):
        for ii, _ in enumerate(supercell):
            if ii not in matches:
                print(f"Missing: {ii} {supercell.positions[ii]}")

    assert len(np.unique(indices, axis=0)) == len(supercell), (indices, len(supercell))

    timer(f"matched {len(matches)} positions in supercell and primitive cell")

    return indices


def get_lattice_points(atoms, supercell, fortran=True, **kwargs):
    """ wrap _get_lattice_points """
    return _get_lattice_points(atoms.cell, supercell.cell, fortran=fortran, **kwargs)


def get_supercell_positions(atoms, supercell):
    """ find all positions in the supercell including the boundary atom """

    timer = Timer()
    tol = tolerance

    lattice_points, lattice_points_ext = get_lattice_points(atoms, supercell)


def _get_lattice_points(
    lattice, superlattice, tolerance=1e-5, sort=True, fortran=True, verbose=False
):
    """
        S = M . L

        M = supercell_matrix """

    timer = Timer()
    tol = tolerance

    inv_lattice = la.inv(lattice)
    inv_superlattice = la.inv(superlattice)

    supercell_matrix = np.round(superlattice @ inv_lattice).astype(int)

    # How many lattice points are to be expected?
    n_lattice_points = int(np.round(np.linalg.det(supercell_matrix)))

    # Maximum number of iterations:
    max_iterations = abs(supercell_matrix).sum()

    if verbose:
        print(f"Maximum number of iterations: {max_iterations}")
        print(f"\nSupercell matrix:             \n{supercell_matrix}")
        print(f"\nlattice:                      \n{lattice}")
        print(f"\ninv_lattice:                  \n{inv_lattice}\n")

    # maximal distance = diagonal of the cell
    # points generated beyond this won't lie inside the supercell
    dmax = 2.5 * np.linalg.norm(superlattice.sum(axis=1))

    if fortran:
        all_lattice_points = sc.supercell.find_lattice_points(
            lattice.T, inv_superlattice.T, n_lattice_points, max_iterations, tolerance
        ).T
        lattice_points = all_lattice_points[:n_lattice_points]
        lattice_points_extended = [
            p
            for p in all_lattice_points[n_lattice_points:]
            if sum(p) > -30000 + tolerance
        ]

    else:
        # find lattice points by enumeration
        lattice_points = []
        lattice_points_extended = []

        for (n1, n2, n3) in product(
            range(-max_iterations, max_iterations + 1), repeat=3
        ):

            lp = [n1, n2, n3] @ lattice

            if la.norm(lp) > dmax:
                continue

            frac_lp = fractional(lp, superlattice)

            # Check if is inside supercell [0, 1) and discard if no
            if (np.all(np.array(frac_lp) > -tolerance)) and (
                np.all(np.array(frac_lp) < 1 - tolerance)
            ):
                lattice_points.append(lp)

            # Check if is inside extended supercell [0, 1] and discard if no
            elif (np.all(np.array(frac_lp) > -tolerance)) and (
                np.all(np.array(frac_lp) < 1 + tolerance)
            ):
                lattice_points_extended.append(lp)

    assert len(np.unique(lattice_points, axis=0)) == n_lattice_points, (
        len(np.unique(lattice_points, axis=0)),
        n_lattice_points,
        lattice_points[:3],
    )

    # assert len(np.unique(lattice_points_extended, axis=0)) == 8 * n_lattice_points, (
    #     len(np.unique(lattice_points_extended, axis=0)),
    #     8 * n_lattice_points,
    #     lattice_points_extended[:3],
    # )

    timer(
        f"found {len(lattice_points)} ({len(lattice_points_extended)}) lattice points"
    )

    if sort:
        lattice_points = np.asarray(sort_lattice_points(lattice_points))
    else:
        lattice_points = np.asarray(lattice_points)

    # find multiplicities of the extended lattice points
    lattice_points_ext_w_multiplicites = []
    for lp in lattice_points:

        frac_lp = fractional(lp, superlattice)

        elp_mult = []

        for elp in lattice_points_extended:
            frac_elp = fractional(elp, superlattice)

            if la.norm((frac_elp - frac_lp + tol) % 1 % 1 - tol) < tol:
                elp_mult.append(elp)

        lattice_points_ext_w_multiplicites.append(elp_mult)

    # len_elps = [len(s) for s in lattice_points_ext_w_multiplicites]
    # assert sum(len_elps) == 8 * n_lattice_points, len_elps

    return lattice_points, lattice_points_ext_w_multiplicites


def sort_lattice_points(lattice_points, tol=1e-5):
    """ sort according to x, y, z coordinates and finally length """

    return sorted(lattice_points, key=lambda x: la.norm(x + [0, 2 * tol, 4 * tol]))


def fractional(positions, lattice):
    """ compute fractioal components in terms of lattice

            r = r_frac . lattice

        =>  r_frac = r . inv_lattice

        """

    frac = positions @ la.inv(lattice)

    return frac
