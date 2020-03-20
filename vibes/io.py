from pathlib import Path

from ase.io import read as ase_read

from vibes.spglib.wrapper import get_symmetry_dataset
from vibes.structure.io import inform  # noqa: F401


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


def parse_force_constants(fc_file, **kwargs):
    """parse either phonopy FORCE_CONSTANTS or tdep infile.forceconstants"""
    file = Path(fc_file)

    name = str(file).lower()

    if "force_constants" in name or "fc2" in name:
        from vibes.phonopy.utils import parse_phonopy_force_constants

        return parse_phonopy_force_constants(file, **kwargs)
    elif name.endswith(".forceconstant"):
        from vibes.tdep.wrapper import parse_tdep_forceconstant

        return parse_tdep_forceconstant(file, **kwargs)

    elif name.endswith(".forceconstant_remapped"):
        from vibes.tdep.wrapper import parse_tdep_remapped_forceconstant

        return parse_tdep_remapped_forceconstant(file, **kwargs)

    elif ".dat" in name or name == "FORCE_CONSTANTS_remapped":
        import numpy as np

        return np.loadtxt(fc_file)

    else:
        raise RuntimeError(f"{file} is neither phonopy nor tdep force constants")
