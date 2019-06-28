from ase.atoms import Atoms
import numpy as np


def to_phonopy_atoms(structure, wrap=False):
    """Convert ase.atoms.Atoms to PhonopyAtoms

    Parameters
    ----------
    structure: ase.atoms.Atoms
        Atoms to convert
    wrap: bool
        If True wrap the scaled positions

    Returns
    -------
    phonopy_atoms: PhonopyAtoms
        The PhonopyAtoms for the same structure as atoms
    """
    from phonopy.structure.atoms import PhonopyAtoms

    phonopy_atoms = PhonopyAtoms(
        symbols=structure.get_chemical_symbols(),
        cell=structure.get_cell(),
        masses=structure.get_masses(),
        positions=structure.get_positions(wrap=wrap),
    )
    return phonopy_atoms


def to_spglib_cell(structure):
    """Convert ase.atoms.Atoms to spglib cell

    Parameters
    ----------
    structure: ase.atoms.Atoms
        Atoms to convert

    Returns
    -------
    lattice: np.ndarray
        The lattice vectors of the cell
    positions: np.ndarray
        The scaled positions of the cell
    number: np.ndarray
        The atomic number of all atoms in the cell
    """
    lattice = structure.cell
    positions = structure.get_scaled_positions()
    number = structure.get_atomic_numbers()
    return (lattice, positions, number)


def to_Atoms_db(structure, info=None, pbc=True):
    """Convert structure to ase.atoms.Atoms without masses, and more accurate positions/lattice vectors

    Parameters
    ----------
    structure: PhonopyAtoms
        The structure to convert
    info: dict
        Additional information to include in atoms.info
    pbc: bool
        True if the structure is periodic

    Returns
    -------
    atoms: ase.atoms.Atoms
        The ASE representation of the material
    """

    if info is None:
        info = {}

    if structure is None:
        return None

    atoms_dict = {
        "symbols": structure.get_chemical_symbols(),
        "cell": np.round(structure.get_cell(), 14),
        "positions": np.round(structure.get_positions(), 14),
        "pbc": pbc,
        "info": info,
    }

    atoms = Atoms(**atoms_dict)

    return atoms


def to_Atoms(structure, info=None, pbc=True):
    """Convert structure to ase.atoms.Atoms

    Parameters
    ----------
    structure: PhonopyAtoms
        The structure to convert
    info: dict
        Additional information to include in atoms.info
    pbc: bool
        True if the structure is periodic

    Returns
    -------
    atoms: ase.atoms.Atoms
        The ASE representation of the material
    """

    if info is None:
        info = {}

    if structure is None:
        return None

    atoms_dict = {
        "symbols": structure.get_chemical_symbols(),
        "cell": structure.get_cell(),
        "masses": structure.get_masses(),
        "positions": structure.get_positions(),
        "pbc": pbc,
        "info": info,
    }

    atoms = Atoms(**atoms_dict)

    return atoms
