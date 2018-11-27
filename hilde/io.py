from ase.io import read as ase_read
from hilde.spglib.wrapper import get_symmetry_dataset
from hilde.structure.io import inform


def get_info_str(atoms, spacegroup=None):
    """ encode atoms.info as string """

    info_strings = [f"Number of atoms:     {len(atoms)}"]

    if spacegroup is not None:
        info_strings += [f"Spacegroup:          {spacegroup}"]

    for key, val in atoms.info.items():
        info_strings.append(f"{key + ':':20s} {val}")

    return info_strings


def read(fname, format="aims"):
    """ wrap ase.io.read """

    atoms = ase_read(fname, format=format)

    return atoms


def write(atoms, fname, format="aims", spacegroup=None, **kwargs):
    """ wrap ase.io.write """

    if format == "aims":
        atoms.write(
            fname, info_str=get_info_str(atoms, spacegroup), format=format, **kwargs
        )

    return True
