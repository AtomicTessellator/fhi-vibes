""" a light wrapper for spglib """
import collections

import numpy as np
import spglib
from ase.atoms import Atoms

from vibes.structure.convert import to_spglib_cell


def cell_to_Atoms(
    lattice: list, scaled_positions: list, numbers: list, info: dict = None
) -> Atoms:
    """convert from spglib cell to Atoms

    Args:
        lattice: The lattice vectors of the structure
        scaled_positions: The scaled positions of the atoms
        numbers: Atomic numbers of all the atoms in the cell
        info: additional information on the structure

    Returns:
        The ASE Atoms object representation of the material

    """
    atoms_dict = {
        "cell": lattice,
        "scaled_positions": scaled_positions,
        "numbers": numbers,
        "pbc": True,
        "info": info,
    }

    return Atoms(**atoms_dict)


def get_symmetry_dataset(
    atoms: Atoms, index_maps: bool = False, symprec: float = 1e-5
) -> collections.namedtuple:
    """return the spglib symmetry dataset

    Args:
        atoms: The structure to get the dataset of (usually primitive cell)
        index_maps: optionally include index maps for the symmetry operations (expensive)

    Returns:
        The symmetry_dataset for the structure including cartesian rotations and
        translations as namedtuple

    """
    spg_cell = to_spglib_cell(atoms)
    data = spglib.get_symmetry_dataset(spg_cell, symprec=symprec)

    uwcks, count = np.unique(data["wyckoffs"], return_counts=True)
    data["wyckoffs_unique"] = [(w, c) for (w, c) in zip(uwcks, count)]

    # add euivalent atoms
    ats, count = np.unique(data["equivalent_atoms"], return_counts=True)
    data["equivalent_atoms_unique"] = [(a, c) for (a, c) in zip(ats, count)]

    # add cartesian rotations and translations
    cell = atoms.cell
    rotations = data["rotations"]
    rotations_cart = [cell.T @ (rot @ cell.reciprocal()) for rot in rotations]
    rotations_cart = np.array(rotations_cart).round(decimals=15)  # remove 0s
    data["rotations_cartesian"] = rotations_cart

    translations = data["translations"]
    translations_cart = [trl @ cell for trl in translations]
    translations_cart = np.array(translations_cart).round(decimals=15)  # remove 0s
    data["translations_cartesian"] = translations_cart

    if index_maps:  # attach index map
        kw = {
            "rotations_cartesian": rotations_cart,
            "translations_cartesian": translations_cart,
        }
        data["index_maps"] = _get_index_maps(atoms, **kw)

    return collections.namedtuple("spg_dataset", data.keys())(**data)


def map_unique_to_atoms(atoms: Atoms, symprec: float = 1e-5) -> list:
    """map each symmetry unique atom to other atoms as used by phonopy PDOS

    Args:
        atoms: The structure to get the dataset of

    Returns:
        The mapping of symmetry unique atoms to other atoms

    """
    ds = get_symmetry_dataset(atoms, symprec=symprec)

    uniques = np.unique(ds.equivalent_atoms)

    mapping = [[] for _ in range(len(uniques))]

    for ii, index in enumerate(ds.equivalent_atoms):
        for jj, unique in enumerate(uniques):
            if index == unique:
                mapping[jj].append(ii)

    return mapping


def get_spacegroup(atoms: Atoms, symprec: float = 1e-5) -> str:
    """return spglib spacegroup

    Args:
        atoms: The structure to get the dataset of

    Returns:
        The spglib space group

    """
    return spglib.get_spacegroup(to_spglib_cell(atoms), symprec=symprec)


def refine_cell(atoms: Atoms, symprec: float = 1e-5) -> Atoms:
    """refine the structure

    Args:
        atoms: The structure to get the dataset of

    Returns:
        The refined structure of atoms

    """
    lattice, scaled_positions, numbers = spglib.refine_cell(
        to_spglib_cell(atoms), symprec
    )

    return cell_to_Atoms(lattice, scaled_positions, numbers)


def standardize_cell(
    atoms: Atoms,
    to_primitive: bool = False,
    no_idealize: bool = False,
    symprec: float = 1e-5,
) -> Atoms:
    """wrap spglib.standardize_cell

    Args:
        atoms: The structure to get the dataset of
        to_primitive: If True go to the primitive cell
        no_idealize: If True do not idealize the cell

    Returns:
        The standardized structure of atoms

    """
    cell = to_spglib_cell(atoms)
    args = spglib.standardize_cell(cell, to_primitive, no_idealize, symprec)

    return cell_to_Atoms(*args)


def _get_index_maps(
    atoms: Atoms, rotations_cartesian, translations_cartesian
) -> np.ndarray:
    """create index maps for structure in atoms for the given symmetry operations

    Args:
        atoms: the structure, e.g., primitive cell
        rotations: rotations matrices IN CARTESIAN COORDS
        translations: translation vectors IN CARTESIAN COORDS

    Returns:
        array containing index maps for each symmetry operation (rot, trl)

    """
    from vibes.helpers.supercell import map_indices

    pos_cart = atoms.positions.copy()

    index_maps = np.zeros([len(rotations_cartesian), len(atoms)], dtype=int)

    for ii, (rc, tc) in enumerate(zip(rotations_cartesian, translations_cartesian)):
        atoms_transformed = atoms.copy()
        atoms_transformed.positions = (rc @ pos_cart.T).T + tc

        index_map = map_indices(atoms, atoms_transformed)
        index_maps[ii, :] = index_map

    return index_maps


def get_index_maps(atoms: Atoms, spg_data: collections.namedtuple) -> np.ndarray:
    """return index maps for structure in Atoms and symmetry operations in spg_data

    Args:
        atoms: the structure
        spg_data: spg data obtained from `get_symmetry_dataset`

    Returns:
        the index maps

    """
    kw = {
        "rotations_cartesian": spg_data.rotations_cartesian,
        "translations_cartesian": spg_data.translations_cartesian,
    }
    return _get_index_maps(atoms, **kw)


def get_ir_reciprocal_mesh(
    mesh: np.ndarray,
    atoms: Atoms,
    monkhorst: bool = True,
    symprec: float = 1e-5,
    eps: float = 1e-9,
) -> collections.namedtuple:
    """wrapper for spglib.get_ir_reciprocal_mesh

    Args:
        mesh: array specifying number of points per axis
        atoms: structure that determines symmetry
        monkhorst: return monkhorst-pack-style grid (gamma incl. when odd grid number)

    Returns:
        (points, weights): frac. grid points in interval [-0.5, 0.5), their weights

    """
    mesh = np.asarray(mesh)
    spg_cell = to_spglib_cell(atoms)

    if monkhorst:  # apply shift
        is_shift = -(np.asarray(mesh) % 2) + 1
    else:
        is_shift = np.zeros(3)

    mapping, grid = spglib.get_ir_reciprocal_mesh(
        mesh, cell=spg_cell, is_shift=is_shift
    )

    map2ir, weights = np.unique(mapping, return_counts=True)
    ir_points = grid[map2ir].astype(float) / mesh  # to frac. coords
    ir_points += 0.5 * is_shift / mesh  # restore shift

    # map to [-0.5, 0.5)
    ir_points = ((ir_points + 0.5 + eps) % 1 - 0.5 - eps).round(decimals=14)
    ir_points_cart = atoms.cell.reciprocal().cartesian_positions(ir_points)

    data = {"points": ir_points, "points_cartesian": ir_points_cart, "weights": weights}

    IrReciprocalMesh = collections.namedtuple("ir_reciprocal_mesh", data.keys())

    return IrReciprocalMesh(**data)
