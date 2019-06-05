""" tools for converting atoms objects to json representations """


import json
from pathlib import Path
import numpy as np
from ase.db.row import atoms2dict as ase_atoms2dict
from ase.io.jsonio import MyEncoder
from ase.calculators.calculator import all_properties
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import voigt_6_to_full_3x3_stress


def atoms2dict(atoms):
    """ Converts an Atoms object into a dict """

    atoms_dict = {
        "symbols": [f"{sym}" for sym in atoms.symbols],
        "masses": atoms.get_masses().tolist(),
    }

    # if periodic system, append lattice (before positions)
    if any(atoms.pbc):
        atoms_dict.update({"cell": atoms.cell.tolist()})

    atoms_dict.update({"positions": atoms.positions.tolist()})

    if atoms.get_velocities() is not None:
        atoms_dict.update({"velocities": atoms.get_velocities().tolist()})

    if atoms.info != {}:
        atoms_dict.update({"info": atoms.info})

    return atoms_dict


def calc2dict(calc):
    """ Converts an ase calculator calc into a dict"""

    if calc is None:
        return {}
    elif isinstance(calc, dict):
        return calc

    params = calc.todict()
    for key, val in params.items():
        if isinstance(val, tuple):
            params[key] = list(val)

    calc_dict = {"calculator": calc.__class__.__name__, "calculator_parameters": params}
    if hasattr(calc_dict, "command"):
        calc_dict.update({"command": calc.command})

    return calc_dict


def input2dict(atoms, calc=None, primitive=None, supercell=None, settings=False):
    """ convert metadata information to plain dict

    Returns:
        {'calculator': calc, 'atoms': atoms} """

    # structure
    atoms_dict = atoms2dict(atoms)

    # calculator
    if calc is None:
        calc = atoms.calc

    calc_dict = calc2dict(calc)

    input_dict = {"calculator": calc_dict, "atoms": atoms_dict}

    if primitive:
        input_dict.update({"primitive": atoms2dict(primitive)})

    if supercell:
        input_dict.update({"supercell": atoms2dict(supercell)})

    # save the configuration
    if settings:
        from hilde.settings import Settings

        settings_dict = dict(Settings())

        input_dict.update({"settings": settings_dict})

    return input_dict


def results2dict(atoms, calc=None, append_cell=False):
    """ extract information from atoms and calculator and convert to plain dict """

    if calc is None:
        calc = atoms.calc

    if atoms.info:
        atoms_dict = {"info": atoms.info}
    else:
        atoms_dict = {}

    # if periodic system, append lattice
    if append_cell and any(atoms.pbc):
        atoms_dict.update({"cell": atoms.cell.tolist()})

    # add positions
    atoms_dict.update({"positions": atoms.positions.tolist()})

    if atoms.get_velocities() is not None:
        atoms_dict.update({"velocities": atoms.get_velocities().tolist()})

    # calculated values
    calc_dict = {}

    # convert stress to 3x3 if present
    if "stress" in calc.results:
        stress = calc.results["stress"]
        if len(stress) == 6:
            calc.results["stress"] = voigt_6_to_full_3x3_stress(stress)

    # convert numpy arrays into ordinary lists
    for key, val in calc.results.items():
        if isinstance(val, np.ndarray):
            calc_dict[key] = val.tolist()
        elif isinstance(val, np.float):
            calc_dict[key] = float(val)
        else:
            calc_dict[key] = val

    return {"atoms": atoms_dict, "calculator": calc_dict}


def dict2atoms(atoms_dict, calc_dict=None):
    """ convert dictionaries into atoms and calculator objects """

    pbc = False
    if "cell" in atoms_dict:
        pbc = True

    try:
        velocities = atoms_dict.pop("velocities")
    except KeyError:
        velocities = None
    atoms = Atoms(**atoms_dict, pbc=pbc)

    if velocities is not None:
        atoms.set_velocities(velocities)

    # Calculator
    if calc_dict is not None:
        results = {}
        if "results" in calc_dict:
            results = calc_dict.pop("results")

        calc = SinglePointCalculator(atoms, **results)
        if "calculator" in calc_dict:
            calc.name = calc_dict["calculator"].lower()
        if "calculator_parameters" in calc_dict:
            calc.parameters.update(calc_dict["calculator_parameters"])
        if "command" in calc_dict:
            calc.command = calc_dict["command"]
    else:
        calc = None

    atoms.calc = calc

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
    atomsdict = ase_atoms2dict(atoms)

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
