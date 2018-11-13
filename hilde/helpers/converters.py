""" tools for conerting atoms objects to json representations """


import json
from pathlib import Path
import numpy as np
from ase.db.row import atoms2dict, AtomsRow
from ase.io.jsonio import MyEncoder
from ase.calculators.calculator import all_properties

from hilde.structure import pAtoms


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
    try:
        atoms = pAtoms(AtomsRow(atoms_dict).toatoms(attach_calculator=True))
    except AttributeError:
        atoms = pAtoms(AtomsRow(atoms_dict).toatoms(attach_calculator=False))

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


def get_json(obj):
    "Return json representation of obj"
    return json.dumps(obj, cls=MyEncoder, sort_keys=True)


def atoms2json(
    atoms, ignore_results=False, ignore_keys=["unique_id"], ignore_calc_params=[]
):
    """ return json representation of atoms and calculator objects.
        possibility to remove certain keys from the atoms dictionary, e.g. for hashing
    """

    # dictionary contains all the information in atoms object
    atomsdict = atoms2dict(atoms)

    # remove unwanted keys from atomsdict
    for name in ignore_keys:
        if name in atomsdict:
            atomsdict.pop(name)

    calcdict = {}

    # move physical properties from atomsdict to calcdict if they are wanted
    for key in all_properties:
        if key in atomsdict:
            value = atomsdict.pop(key)
            if not ignore_results:
                calcdict[key] = value

    # clean calculator entries
    if "calculator_parameters" in atomsdict:
        calculator_params = atomsdict["calculator_parameters"]
        for name in [key for key in calculator_params if key in ignore_calc_params]:
            calculator_params.pop(name)

        if "species_dir" in calculator_params:
            calculator_params["species_dir"] = Path(
                calculator_params["species_dir"]
            ).parts[-1]

    for name in ["calculator", "calculator_parameters"]:
        if name in atomsdict:
            calcdict[name] = atomsdict.pop(name)

    return get_json(atomsdict), get_json(calcdict)
