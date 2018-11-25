"""
tools for converting atoms and calculator objects to and from different representations
"""
import numpy as np
from ase.db.row import atoms2dict, AtomsRow

from hilde.structure import pAtoms
from hilde.calculators.aims_calc import Aims


def calc2dict(calc):
    """ Converts an ase calculator calc into a dict"""
    if calc is None:
        return {}
    calc_dict = {}
    calc_dict["calculator"] = calc.name.lower()
    calc_dict["calculator_parameters"] = calc.todict()
    try:
        calc_dict["command"] = calc.command
    except:
        pass
    calc_dict["results"] = calc.results
    return calc_dict


def patoms2dict(atoms):
    """
    Converts a pAtoms object into a dict
    Args:
        atoms: pAtoms or Atoms object
            The pAtoms or Atoms object to be converted into a dictionary
    Returns: atoms_dict (dict)
        The dictionary of atoms
    """
    if atoms is None:
        return atoms
    atoms_dict = atoms2dict(atoms)

    # add information that is missing after using ase.atoms2dict
    atoms_dict["info"] = atoms.info

    if hasattr(atoms, "symmetry_block"):
        atoms_dict["sym_block"] = atoms.symmetry_block

    # attach calculator
    for key, val in calc2dict(atoms.calc).items():
        atoms_dict[key] = val

    return atoms_dict


def dict2patoms(atoms_dict):
    """
    Converts a dict into a pAtoms object
    Args:
        atoms_dict: dict
            A dictionary representing the pAtoms object
    Returns: pAtoms
        The corresponding pAtoms object
    """
    if atoms_dict is None:
        return None
    try:
        atoms = pAtoms(AtomsRow(atoms_dict).toatoms(attach_calculator=True))
    except AttributeError:
        atoms = pAtoms(AtomsRow(atoms_dict).toatoms(attach_calculator=False))

    if "calculator" in atoms_dict and atoms_dict["calculator"] == "aims":
        atoms.calc = Aims(
            aims_command=atoms_dict["command"], **atoms_dict["calculator_parameters"]
        )

    # Attach missing information
    if "info" in atoms_dict:
        atoms.info = atoms_dict["info"]
    if "sym_block" in atoms_dict:
        atoms.symmetry_block = atoms_dict["sym_block"]
    if "command" in atoms_dict:
        atoms.calc.command = atoms_dict["command"]
    if "results" in atoms_dict:
        atoms.calc.results = atoms_dict["results"]

    # attach calculator
    if atoms.calc:
        for key, val in atoms.calc.results.items():
            if isinstance(val, list):
                atoms.calc.results[key] = np.array(val)

    return atoms
