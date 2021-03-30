""" helps to find lattice points in supercell, match positions to images in the
unit cell etc. """

import collections
from itertools import product

import numpy as np
import scipy.linalg as la
from ase.cell import Cell

from vibes.ase.calculators.fc import get_smallest_vectors  # noqa: F401
from vibes.helpers.supercell import supercell as sc
from vibes.helpers.utils import Timer


_prefix = "lattice_points"


def get_lattice_points(
    cell: np.ndarray,
    supercell: np.ndarray,
    extended: bool = True,
    tolerance: float = 1e-5,
    sort: bool = True,
    fortran: bool = True,
    decimals: int = 16,
    verbose: bool = False,
) -> (np.ndarray, list):
    """
        S = M . L

        M = supercell_matrix

    Args:
        cell: lattice matrix of primitive cell
        supercell: lattice matrix of supercell
        extended: return extendes lattice points at supercell boundary
        tolerance: tolerance used to detect multiplicities
        sort: If True sort results
        fortran: If True use the Fortran routines
        decimals: how many digits to round to
        verbose: If True print more information on the console

    Returns:
        lattice_points: Array of lattice points in the supercell
        lattice_points_ext_w_multiplicites: lattice points including multiplicities
    """

    timer = Timer(prefix=_prefix, verbose=verbose)
    tol = tolerance

    # check if cells are arrays
    pcell = Cell(cell)
    scell = Cell(supercell)

    lattice = pcell[:]
    superlattice = scell[:]

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

            frac_lp = scell.scaled_positions(lp)

            # Check if is inside supercell [-0.5, 0.5) and discard if no
            if (np.all(np.array(frac_lp) > -0.5 - tolerance)) and (
                np.all(np.array(frac_lp) < 0.5 - tolerance)
            ):
                lattice_points.append(lp)

            # Check if is inside extended supercell [-0.5, 0.5] and discard if no
            elif (np.all(np.array(frac_lp) > -0.5 - tolerance)) and (
                np.all(np.array(frac_lp) < 0.5 + tolerance)
            ):
                lattice_points_extended.append(lp)

    assert len(np.unique(lattice_points, axis=0)) == n_lattice_points, (
        len(np.unique(lattice_points, axis=0)),
        n_lattice_points,
        lattice_points[:3],
    )

    nlp = len(lattice_points)
    nlpe = len(lattice_points_extended)
    timer(f"found {nlp} ({nlpe}) lattice points")

    if sort:
        lattice_points = sort_lattice_points(lattice_points)

    lattice_points = np.around(lattice_points, decimals=decimals)

    if not extended:
        return lattice_points

    # find multiplicities of the extended lattice points
    lattice_points_ext_w_multiplicites = []
    for lp in lattice_points:

        frac_lp = scell.scaled_positions(lp)

        elp_mult = [lp]

        for elp in lattice_points_extended:
            frac_elp = scell.scaled_positions(elp)

            if la.norm((frac_elp - frac_lp + tol) % 1 % 1 - tol) < tol:
                elp_mult.append(elp)

        elp_mult = np.around(elp_mult, decimals=decimals)
        lattice_points_ext_w_multiplicites.append(elp_mult)

    return lattice_points, lattice_points_ext_w_multiplicites


def sort_lattice_points(lattice_points, tol=1e-5):
    """sort according to x, y, z coordinates and finally length

    Parameters
    ----------
    lattice_points: np.ndarray
        The list of lattice points in the supercell
    tol: float
        tolerance for small numbers

    Returns
    -------
    np.ndarray
        sorted lattice point list
    """

    return sorted(lattice_points, key=lambda x: la.norm(x + [0, 2 * tol, 4 * tol]))


def map_L_to_i(indeces):
    """ Map to atoms belonging to specific lattice point

    Parameters
    ----------
    indeces: list
        map from u_I in supercell to u_iL w.r.t to primitive cell and lattice point

    Returns
    -------
    np.ndarray
        list of masks that single out the atoms that belong to specific lattice point
    """

    n_lattice_points = max(i[1] for i in indeces) + 1
    mappings = []
    for LL in range(n_lattice_points):
        mappings.append([idx[1] == LL for idx in indeces])
    return mappings


def map_I_to_iL(
    in_atoms, in_supercell, lattice_points=None, tol=1e-5, verbose=True
) -> tuple:
    """Map from supercell index I to (i, L), i is the unit cell index and L lattice p.

    Args:
        in_atoms (ase.atoms.Atoms): primitive cell
        in_supercell (ase.atoms.Atoms): supercell
        lattice_points (list, optional): list of pre-computed lattice points L
        tol (float, optional): tolerance for wrapping

    Returns:
        list, list: I_to_iL map, inverse map

    """
    timer = Timer(prefix=_prefix, verbose=verbose)

    atoms = in_atoms.copy()
    supercell = in_supercell.copy()
    atoms.wrap()
    supercell.wrap()

    if lattice_points is None:
        lattice_points, _ = get_lattice_points(atoms.cell, supercell.cell)

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
                if la.norm(spos - pos - lp) < tol:
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
                fpos = supercell.cell.scaled_positions(spos - pos - lp)
                if la.norm((fpos + tol) % 1 - tol) < tol:
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

    inv = _map_iL_to_I(indices)

    return np.array(indices), np.array(inv)


def _map_iL_to_I(I_to_iL_map):
    """ map (i, L) back to supercell index I

    Parameters
    ----------
    I_to_iL_map: np.ndarray
        Map from I to iL

    Returns
    -------
    np.ndarray
        Map back from primitive cell index/lattice point to supercell index

    Raises
    ------
    AssertionError
        If iL2I[I2iL[II][0], I2iL[II][1]] does not equal II
    """

    I2iL = np.array(I_to_iL_map)

    n_atoms = int(I2iL[:, 0].max() + 1)
    n_lps = int(I2iL[:, 1].max() + 1)

    iL2I = np.zeros([n_atoms, n_lps], dtype=int)

    for II, (ii, LL) in enumerate(I_to_iL_map):
        iL2I[ii, LL] = II

    # sanity check:
    for II in range(n_atoms * n_lps):
        iL = I2iL[II]
        I = iL2I[iL[0], iL[1]]
        assert II == I, (II, iL, I)

    return iL2I.squeeze()


def get_commensurate_q_points(
    cell: np.ndarray,
    supercell: np.ndarray,
    fractional: bool = False,
    decimals: int = 14,
    tolerance: float = 1e-5,
) -> np.ndarray:
    """For a commensurate q_points we have

        exp( 2*pi q . L_k ) = 1 for any k and L_k being the supercell lattice vectors

        in other workds, q is a linear combination of G_k, where G_k are the inverse
        lattice vectors of the supercell lattice. Only those are counted which fit into
        the inverse lattice of the primitive cell.
        This means we have to call lattice_points.get_lattice_points with the inverse
        lattices.

    Args:
        cell: cell matrix of primitive cell
        supercell: cell matrix of supercell
        fractional: return in fractional coords
        decimals: use to clean zeros
        tolerance: tolearnce used to detect multiplicities

    Returns:
        List of commensurate q_points

    """
    lattice = cell
    superlattice = supercell

    inv_lattice = la.inv(lattice).T
    inv_superlattice = la.inv(superlattice).T

    inv_lattice_points = get_lattice_points(
        inv_superlattice, inv_lattice, extended=False, tolerance=tolerance
    )

    if fractional:
        filps = la.solve(inv_lattice.T, np.transpose(inv_lattice_points)).T
        assert np.allclose(filps, inv_lattice_points @ cell.T), "FIX ME"
        inv_lattice_points = filps.round(decimals=decimals)

    return inv_lattice_points


def get_unit_grid_extended(
    q_points_frac: np.ndarray, only_gamma: bool = False, tol: float = 1e-9
) -> tuple:
    """map q-points to units cube [0, 1] adding boundary elements for interpolation

    Args:
        q_points_frac: q-points to be extended, e.g. in interval [-0.5, 0.5)
        only_gamma: only include gamma point (for symmetry-reduced grids).
                    REM: assume gamma point is the first point

    Returns:
        (grid, grid_extended, map_to_extended):
            grid on unit cube [0, 1), extended to [0, 1], map from extended to orig.
    """
    # map to [0, 1]
    q_points_frac_unit = (q_points_frac + tol) % 1 - tol

    map_to_extended = []
    extended_q_points_unit = []
    for iq, q in enumerate(q_points_frac_unit):
        if (abs(q) < tol).any():
            # print(iq, q)
            # create all reciprocal lattice vector displacements from 0 -> 1
            mask = abs(q) < tol  # find elements that are 0
            rge = np.ones(3)
            for ii in range(3):  # create range for ndindex
                if mask[ii]:
                    rge[ii] = 2

            for ii, jj, kk in np.ndindex(*rge.astype(int)):
                if ii == jj == kk == 0:
                    continue
                disp = np.array([ii, jj, kk])
                new_q = q + disp
                # print(new_q)
                extended_q_points_unit.append(new_q)
                map_to_extended.append(iq)

        if only_gamma:
            assert la.norm(q) < tol, "first point was not gamma, not implemented, ADD!"
            break

    q_points_frac_unit_extended = np.concatenate(
        (q_points_frac_unit, extended_q_points_unit)
    )
    map_to_extended = np.concatenate((range(len(q_points_frac)), map_to_extended))
    # w_sq_extended = np.concatenate((dmx.w_sq, np.asarray(extended_w_sq).T), axis=1)

    return collections.namedtuple(
        "unit_grid_extended", ("points", "points_extended", "map2extended"),
    )(q_points_frac_unit, q_points_frac_unit_extended, map_to_extended)
