""" helpers for working with supercell """

from itertools import product
import numpy as np
import scipy.linalg as la

from ase import Atoms
from ase.spacegroup import get_spacegroup
from hilde.structure.misc import get_sysname
from hilde.helpers.utils import Timer
from hilde.helpers.geometry import get_cubicness
from hilde.helpers.lattice import fractional
from hilde.helpers.numerics import clean_matrix
from hilde.helpers.warnings import warn
from hilde.phonopy.wrapper import preprocess
from . import supercell as sc


def get_lattice_points(
    cell, supercell, tolerance=1e-5, sort=True, fortran=True, verbose=False
):
    """
        S = M . L

        M = supercell_matrix

    Args:
        cell (ndarray): lattice matrix of primitive cell
        supercell (ndarray): lattice matrix of supercell
        """

    timer = Timer()
    tol = tolerance

    # check if cell is array
    if isinstance(cell, Atoms):
        warn("DEPRECATED: Please provide 3x3 matrix instead of Atoms object")
        cell = cell.cell
    # check if cell is array
    if isinstance(supercell, Atoms):
        warn("DEPRECATED: Please provide 3x3 matrix instead of Atoms object")
        supercell = supercell.cell

    if verbose:
        print("get_lattice_points()\n--------------------")

    lattice = cell.copy()
    superlattice = supercell.copy()

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


def get_commensurate_q_points(cell, supercell, tolerance=1e-5, **kwargs):
    """ For a commensurate q_points we have

        exp( 2*pi q . L_k ) = 1 for any k and L_k being the supercell lattice vectors

        in other workds, q is a linear combination of G_k, where G_k are the inverse
        lattice vectors of the supercell lattice. Only thos are counted which fit into
        the inverse lattice of the primitive cell.
        This means we have to call lattice_points.get_lattice_points with the inverse
        lattices.

    Args:
        cell (ndarray): cell matrix of primitive cell
        supercell (ndarray): cell matrix of supercell

    """

    # check if cell is array
    if isinstance(cell, Atoms):
        warn("DEPRECATED: Please provide 3x3 matrix instead of Atoms object")
        cell = cell.cell
    # check if cell is array
    if isinstance(supercell, Atoms):
        warn("DEPRECATED: Please provide 3x3 matrix instead of Atoms object")
        supercell = supercell.cell

    lattice = cell.T
    superlattice = supercell.T

    inv_lattice = la.inv(lattice)
    inv_superlattice = la.inv(superlattice)

    inv_lattice_points, _ = get_lattice_points(
        inv_superlattice, inv_lattice, tolerance, **kwargs
    )

    return clean_matrix(inv_lattice_points)


def find_cubic_cell(
    cell, target_size=1, deviation=0.2, lower_limit=-2, upper_limit=2, verbose=False
):
    """ Find a supercell matrix that produces a supercell of given size that is
    as cubic as possible """

    smatrix = sc.supercell.find_optimal_cell(
        cell,
        np.eye(3),
        target_size=target_size,
        deviation=deviation,
        lower_limit=lower_limit,
        upper_limit=upper_limit,
        verbose=verbose,
    )
    return np.asarray(smatrix, dtype=int)


def make_cubic_supercell(atoms, target_size=100, deviation=0.2, limit=2, verbose=False):
    """ Create a supercell of target size that is as cubic as possible.

    Args:
        atoms (Atoms): Input atoms object
        target_size (int): Number of atoms in supercell
        deviation (float): Allowed deviation from target supercell size
        limit (int): limit for expansion about analytic search
        verbose (boolean): be verbose (for debugging)

    Returns:
        (Atoms, np.ndarray): supercell, supercell_matrix

    """

    prim_cell = atoms.copy()

    smatrix = find_cubic_cell(
        cell=prim_cell.cell,
        target_size=target_size / len(prim_cell),
        deviation=deviation,
        lower_limit=-limit,
        upper_limit=limit,
        verbose=verbose,
    )

    supercell = make_supercell(
        prim_cell, smatrix, info={"supercell_matrix": smatrix.flatten().tolist()}
    )

    n_sc = get_spacegroup(supercell).no
    n_at = get_spacegroup(prim_cell).no
    if n_sc != n_at:
        warn("Spacegroup of supercell: " + f"{n_sc} |= {n_at} of reference cell.")

    cub_ness = get_cubicness(supercell.cell)
    if cub_ness < 0.8:
        print(
            "**Warning: Cubicness of supercell is "
            + f"{cub_ness:.3f} ({cub_ness**3:.3f})"
        )
        print(f"**-> Sytems: {get_sysname(prim_cell)}, target size {target_size}")
    return supercell, smatrix


def make_supercell(atoms, supercell_matrix, info={}, tol=1e-5, wrap=True):
    """ Create the lattice points within supercell and attach atoms to each of them
    Args:
        atoms (Atoms): primitive cell as atoms object
        supercell_matrix (ndarray): supercell matrix M with convention A = M . a
        info (dict): attach info dictionary to supercell atoms
        tol (float): numerical tolerance for finding lattice points """

    _, supercell, _ = preprocess(atoms, supercell_matrix)
    supercell.cell = clean_matrix(supercell.cell)
    if wrap:
        supercell.set_scaled_positions(supercell.get_scaled_positions(wrap=True))
    return supercell


def map_indices(atoms1, atoms2, tol=1e-5):
    """ return indices of atoms in atoms1 in atoms2.

    Example:
        atoms1 = [H, O1, O2]
        atoms2 = [O1, H, O2]

        -> map_indices(atoms1, atoms2) = [1, 0, 2]

    For background, see
    https://gitlab.com/flokno/hilde/blob/devel/examples/devel/sort_atoms/sort.ipynb"""

    from hilde.helpers.lattice_points import map_I_to_iL

    _, index_map = map_I_to_iL(atoms2, atoms1, return_inverse=True)

    assert len(np.unique(index_map)) == len(atoms1)

    return index_map
