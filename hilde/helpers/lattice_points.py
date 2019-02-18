""" helps to find lattice points in supercell, match positions to images in the 
unit cell etc. """

from itertools import product
import numpy as np
import scipy.linalg as la
from hilde.helpers.timer import Timer
from hilde.helpers.lattice import fractional
from hilde.helpers.supercell import supercell as sc


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


def map_L_to_i(indeces):
    """ Map to atoms belonging to specific lattice point

        indeces:
            map from u_I in supercell to u_iL w.r.t to primitive cell and lattice point
            i: corresponding atom in primitive cell L: lattice point index
        returns:
            list of masks that single out the atoms in the supercell that belong to
            specific lattice point """

    n_lattice_points = max([i[1] for i in indeces]) + 1
    mappings = []
    for LL in range(n_lattice_points):
        mappings.append([idx[1] == LL for idx in indeces])
    return mappings


def map_I_to_iL(
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

    # prepare the list of indices
    indices = len(supercell) * [(-1, -1)]
    matches = []

    for satom in supercell:
        spos, ssym, jj = satom.position, satom.symbol, satom.index
        for atom in atoms:
            pos, sym, ii = atom.position, atom.symbol, atom.index
            # discard rightaway if not the correct species
            if ssym != sym:
                continue
            for LL, lp in enumerate(lattice_points):
                if la.norm(spos - pos - lp) < tolerance:
                    indices[jj] = (ii, LL)
                    matches.append(jj)
                    break

    # catch possibly unwrapped atoms
    for satom in supercell:
        spos, ssym, jj = satom.position, satom.symbol, satom.index
        if jj in matches:
            continue
        for LL, lp in enumerate(lattice_points):
            for atom in atoms:
                pos, sym, ii = atom.position, atom.symbol, atom.index
                if ssym != sym:
                    continue
                fpos = fractional(spos - pos - lp, supercell.cell)
                tol = tolerance
                if la.norm((fpos + tol) % 1 - tol) < tolerance:
                    indices[jj] = (ii, LL)
                    matches.append(jj)
                    break

    # sanity checks:
    if len(np.unique(matches)) != len(supercell):
        for ii, _ in enumerate(supercell):
            if ii not in matches:
                print(f"Missing: {ii} {supercell.positions[ii]}")

    assert len(np.unique(indices, axis=0)) == len(supercell), (indices, len(supercell))

    # should never arrive here
    assert not any(-1 in l for l in indices), ("Indices found: ", indices)

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

    return lattice_points, lattice_points_ext_w_multiplicites


def sort_lattice_points(lattice_points, tol=1e-5):
    """ sort according to x, y, z coordinates and finally length """

    return sorted(lattice_points, key=lambda x: la.norm(x + [0, 2 * tol, 4 * tol]))
