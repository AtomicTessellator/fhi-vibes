import numpy as np
import scipy.linalg as la
from itertools import product
from hilde.helpers.timer import Timer


def assign_primitive_positions(atoms, supercell, tolerance=1e-5):
    """ write positions in the supercell as (i, L), where i is the respective index
        in the primitive cell and L is the lattice point """

    timer = Timer()

    lattice_points = get_lattice_points(atoms.cell, supercell.cell)

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
        for ii, pos in enumerate(atoms.positions):
            for LL, lp in enumerate(lattice_points):
                for _, LP in enumerate(supercell.cell):
                    if la.norm(spos - pos - lp - LP) < tolerance:
                        indices.append((ii, LL))
                        matches.append(jj)
                        break

    # sanity checks:
    if len(np.unique(matches, axis=0)) != len(supercell):
        for ii, _ in enumerate(supercell):
            if ii not in match:
                print(f"Missing: {ii} {supercell.positions[ii]}")

    assert len(np.unique(indices, axis=0)) == len(supercell), (indices, len(supercell))

    timer(f"matched {len(matches)} positions in supercell and primitive cell")

    return indices


def get_lattice_points(
    lattice, superlattice, tolerance=1e-5, sorted=True, verbose=False
):
    """
        S = M . L

        M = supercell_matrix """

    timer = Timer()

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
        if np.any(np.array(frac_lp) < -tolerance):
            continue

        if np.any(np.array(frac_lp) > 1 - tolerance):
            continue

        # attach to list if passed
        lattice_points.append(lp)

    assert len(np.unique(lattice_points, axis=0)) == n_lattice_points, (
        len(np.unique(lattice_points, axis=0)),
        n_lattice_points,
    )

    timer(f"found {len(lattice_points)} lattice points")

    if sorted:
        return sorted(lattice_points, key=la.norm)
    return lattice_points


def fractional(positions, lattice):
    """ compute fractioal components in terms of lattice

            r = r_frac . lattice

        =>  r_frac = r . inv_lattice

        """

    frac = positions @ la.inv(lattice)

    return frac
