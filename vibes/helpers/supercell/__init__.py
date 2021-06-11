""" helpers for working with supercell """


import numpy as np
from ase import Atoms
from ase.spacegroup import get_spacegroup
from vibes.ase.calculators.fc import get_smallest_vectors  # noqa: F401
from vibes.helpers.geometry import get_cubicness
from vibes.helpers.numerics import clean_matrix
from vibes.helpers.warnings import warn
from vibes.structure.misc import get_sysname


def find_cubic_cell(
    cell, target_size=1, deviation=0.2, lower_limit=-2, upper_limit=2, verbose=False
):
    """Find supercell matrix that produces a cubic-as-possible supercell of given size

    Parameters
    ----------
    cell: np.ndarray
        lattice vectors of the primitive matrix
    target_size: int
        target number of atoms in the supercell
    deviation: float
        acceptable deviation from target size
    lower_limit: int
        lower limit for the elements in the supercell matrix
    upper_limit:int
        upper limit for the elements in the supercell matrix
    verbose: bool
        If True print more information to the console
    """
    from . import supercell as sc

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
    """Create a supercell of target size that is as cubic as possible.

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Input atoms object
    target_size: int
        Number of atoms in supercell
    deviation: float
        Allowed deviation from target supercell size
    limit: int
        limit for expansion about analytic search
    verbose: boolean
        be verbose (for debugging)

    Returns
    -------
    ase.atoms.Atoms
        The supercell
    np.ndarray
        The supercell_matrix

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
        print(f"**-> Systems: {get_sysname(prim_cell)}, target size {target_size}")
    return supercell, smatrix


def make_supercell(atoms, supercell_matrix, info={}, tol=1e-5, wrap=True):
    """Create the lattice points within supercell and attach atoms to each of them

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        primitive cell as atoms object
    supercell_matrix: np.ndarray
        supercell matrix M with convention A = M . a
    info: dict
        attach info dictionary to supercell atoms
    tol: float
        numerical tolerance for finding lattice points

    Returns
    -------
    supercell: ase.atoms.Atoms
        The supercell from atoms and supercell_matrix
    """
    from vibes.phonopy.wrapper import preprocess

    _, supercell, _ = preprocess(atoms, supercell_matrix)
    supercell.cell = clean_matrix(supercell.cell)
    if wrap:
        supercell.set_scaled_positions(supercell.get_scaled_positions(wrap=True))
    return supercell


def map_indices(atoms1: Atoms, atoms2: Atoms, tol: float = 1e-5) -> np.ndarray:
    """ return indices of atoms in atoms1 in atoms2.

    Background Information at
    https://gitlab.com/flokno/vibes/blob/devel/examples/devel/sort_atoms/sort.ipynb

    Example:
        atoms1 = [H, O1, O2]
        atoms2 = [O1, H, O2]

        -> map_indices(atoms1, atoms2) = [1, 0, 2]

    Args:
        atoms1: Structure to get map from
        atoms2: Structure to get map to
        tol: Tolerance to tell if atoms are equivalent

    Returns:
        index_map: Map from atoms1 to atoms2

    """
    from vibes.helpers.lattice_points import map_I_to_iL

    _, index_map = map_I_to_iL(atoms2, atoms1, verbose=False)

    assert len(np.unique(index_map)) == len(atoms1)

    return index_map


def map2prim(primitive: Atoms, supercell: Atoms, tol: float = 1e-5) -> list:
    """Map atoms from supercell to primitive cell and return index map"""
    map2prim = []
    primitive = primitive.copy()
    supercell = supercell.copy()

    # represent new supercell in fractional coords of primitive cell
    supercell_with_prim_cell = supercell.copy()

    supercell_with_prim_cell.cell = primitive.cell.copy()

    primitive.wrap(eps=tol)
    supercell_with_prim_cell.wrap(eps=tol)

    # create list that maps atoms in supercell to atoms in primitive cell
    for a1 in supercell_with_prim_cell:
        diff = primitive.positions - a1.position
        map2prim.append(np.where(np.linalg.norm(diff, axis=1) < tol)[0][0])

    # make sure every atom in primitive was matched equally often
    _, counts = np.unique(map2prim, return_counts=True)
    assert counts.std() == 0, counts

    return map2prim
