""" helpers for working with supercell """

import numpy as np
from ase.spacegroup import get_spacegroup
from ase.build.supercells import make_supercell as ase_make_supercell
from hilde.structure.misc import get_sysname
from hilde.helpers.geometry import get_cubicness
from hilde.helpers.numerics import clean_matrix
from hilde.helpers.warnings import warn
from . import supercell as sc


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


def make_supercell(*args, info={}, **kwargs):
    """ Wrap the make_supercell() function from ase.build """
    supercell = ase_make_supercell(*args, **kwargs)
    supercell.cell = clean_matrix(supercell.cell)
    supercell.info.update(info)

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
