from pathlib import Path

import numpy as np

from ase import Atoms
from ase.formula import Formula
from ase.io import read as ase_read
from ase.spacegroup import get_spacegroup
from vibes.filenames import filenames
from vibes.spglib import get_symmetry_dataset
from vibes.structure.io import inform  # noqa: F401


def get_identifier(atoms: Atoms, fix_names: dict = None) -> dict:
    """Get geometry identifier as dictionary.

    Args:
        atoms: the structure
        fix_names: dict for correcting material names (e.g. OMg -> MgO)
    Returns:
        dict: w/ space_group, n_formula, material (name)

    """
    name = atoms.get_chemical_formula(mode="metal", empirical=True)
    formula = Formula(name)

    count = formula.count()
    nf = sum(count.values())
    sg = get_spacegroup(atoms).no

    if fix_names is not None and name in fix_names:
        name = fix_names[name]

    return {
        "space_group": int(sg),
        "n_formula": int(nf),
        "material": name,
    }


def get_info_str(atoms, spacegroup=False):
    """encode atoms.info as string

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The Atoms object to get the info string for
    spacegroup: bool
        If True add space group information

    Returns
    -------
    info_strings: str
        Teh info string
    """

    info_strings = [f"Number of atoms:     {len(atoms)}"]

    if spacegroup:
        spacegroup = get_symmetry_dataset(atoms).number
        info_strings += [f"Spacegroup:          {spacegroup}"]

    for key, val in atoms.info.items():
        info_strings.append(f"{key + ':':20s} {val}")

    return info_strings


def read(file, format="aims"):
    """wrap ase.io.read

    Parameters
    ----------
    file: str
        The input geometry file
    format: str
        The format of the geometry file

    Returns
    -------
    atoms: ase.atoms.Atoms
        The ASE representation of the structure in file
    """

    atoms = ase_read(file, format=format)

    return atoms


def write(atoms, file, format="aims", spacegroup=False, **kwargs):
    """wrap ase.io.write

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The structure to write to the file
    file: str
        The input geometry file
    format: str
        The format of the geometry file
    spacegroup: bool
        If True add space group information

    Returns
    -------
    bool
        True if completed
    """

    if format == "aims":
        atoms.write(
            file, info_str=get_info_str(atoms, spacegroup), format=format, **kwargs
        )

    return True


def parse_force_constants(fc_file: str, **kwargs) -> np.ndarray:
    """parse either phonopy FORCE_CONSTANTS or numpy array

    Args:
        fc_file: the file with forceconstants

    Returns:
        force_constants as ndarray
    """

    file = Path(fc_file)
    name = file.name

    if ".dat" in name or filenames.fc.phonopy_remapped in name:
        import numpy as np

        return np.loadtxt(fc_file)

    elif filenames.fc.phonopy in name or filenames.fc.phonopy_hdf5 in name:
        from vibes.phonopy.utils import parse_phonopy_force_constants

        return parse_phonopy_force_constants(file, **kwargs)

    else:
        raise RuntimeError(f"{file} type is unkown")
