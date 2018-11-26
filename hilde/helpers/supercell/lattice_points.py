import numpy as np
import scipy.linalg as la
from itertools import product
from hilde.konstanten.numerics import eps, wrap_tol


def assign_primitive_positions(atoms, supercell):
    """ write positions in the supercell as (i, L), where i is the respective index
        in the primitive cell and L is the lattice point """

    lattice_points = sorted(
        find_lattice_points(atoms.cell, supercell.cell), key=la.norm
    )

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
                if la.norm(spos - pos - lp) < wrap_tol:
                    indices.append((ii, LL))
                    matches.append(jj)

    # now check for twice the lattice points, those are probably outliers
    for jj, spos in enumerate(supercell.positions):
        if jj in matches:
            continue
        for ii, pos in enumerate(atoms.positions):
            for LL, lp in enumerate(lattice_points):
                if la.norm(spos - pos - 2 * lp) < wrap_tol:
                    indices.append((ii, 0))
                    matches.append(jj)

    assert len(matches) == len(supercell), (len(matches), len(supercell))

    return indices


def find_lattice_points(lattice, superlattice, verbose=False):
    """
        S = M . L

        M = supercell_matrix """

    inv_lattice = la.inv(lattice)  # .T
    inv_superlattice = la.inv(superlattice)

    supercell_matrix = np.round(superlattice @ inv_lattice).astype(int)

    # How many lattice points are to be expected?
    n_lattice_points = int(np.round(np.linalg.det(supercell_matrix)))

    # Maximum number of iterations:
    max_iterations = abs(supercell_matrix).sum()

    if verbose:
        print(f"Maximum number of iterations: {max_iterations}")

    # maximal distance = diagonal of the cell
    # points generated beyond this won't lie inside the supercell
    dmax = 1.5 * np.linalg.norm(superlattice.sum(axis=1))

    # find lattice points by enumeration
    lattice_points = []
    for (n1, n2, n3) in product(range(-max_iterations, max_iterations + 1), repeat=3):

        lp = [n1, n2, n3] @ lattice

        if la.norm(lp) > dmax:
            continue

        frac_lp = fractional(lp, superlattice)

        # Check if is inside supercell and discard if no
        if np.any(np.array(frac_lp) < -eps):
            continue

        if np.any(np.array(frac_lp) > 1 - eps):
            continue

        # attach to list if passed
        lattice_points.append(lp)

    assert len(lattice_points) == n_lattice_points, (
        len(lattice_points),
        n_lattice_points,
    )

    return lattice_points


def fractional(positions, lattice):
    """ compute fractioal components in terms of lattice

            r = r_frac . lattice

        =>  r_frac = r . inv_lattice

        """

    frac = positions @ la.inv(lattice)

    return frac
