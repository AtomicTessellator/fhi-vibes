""" tools for converting atoms objects to json representations """


import json
from pathlib import Path

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import all_properties, get_calculator_class
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import dict2constraint, voigt_6_to_full_3x3_stress
from ase.db.row import atoms2dict as ase_atoms2dict
from ase.io.jsonio import MyEncoder

from vibes.helpers.lists import expand_list, list_dim, reduce_list
from vibes.konstanten.io import n_yaml_digits


key_symbols = "symbols"
key_masses = "masses"
key_constraints = "constraint"


class NumpyEncoder(json.JSONEncoder):
    """ Decode numerical objects that json cannot parse by default"""

    def default(self, obj):
        """Default JSON encoding

        Parameters
        ----------
        obj
            The object to be encoded

        Returns
        -------
        any
            The JSON friendly version of the object
        """
        if hasattr(obj, "tolist") and callable(obj.tolist):
            return obj.tolist()
        if isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        if isinstance(obj, complex):
            return (float(obj.real), float(obj.imag))
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)


def atoms2dict(atoms, reduce=True, add_constraints=False):
    """Converts an Atoms object into a dict

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The structure to be converted to a dict
    reduce: bool
        use reduced representation of `symbols` and `masses`

    Returns
    -------
    atoms_dict: dict
        The dict representation of atoms
    """

    atoms_dict = {}

    # if periodic system, append lattice (before positions)
    if any(atoms.pbc):
        atoms_dict.update({"pbc": np.asarray(atoms.pbc)})
        atoms_dict.update({"cell": atoms.cell.tolist()})

    atoms_dict.update({"positions": atoms.positions.tolist()})

    if atoms.get_velocities() is not None:
        atoms_dict.update({"velocities": atoms.get_velocities().tolist()})

    # append symbols and masses
    atoms_dict.update(
        {
            key_symbols: reduce_list(atoms.get_chemical_symbols(), reduce=reduce),
            key_masses: reduce_list(atoms.get_masses(), reduce=reduce),
        }
    )

    atoms_dict.update({"info": atoms.info})

    if add_constraints:
        atoms_dict[key_constraints] = [constr.todict() for constr in atoms.constraints]

    return atoms_dict


def calc2dict(calc):
    """Converts an ase calculator calc into a dict

    Parameters
    ----------
    calc: ase.calculators.calulator.Calculator
        The calculator to be converted to a dict

    Returns
    -------
    calc_dict: dict
        The dict representation of calc
    """

    if calc is None:
        return {}
    if isinstance(calc, dict):
        return calc

    params = calc.todict()
    for key, val in params.items():
        if isinstance(val, tuple):
            params[key] = list(val)

    calc_dict = {"calculator": calc.__class__.__name__, "calculator_parameters": params}
    if hasattr(calc_dict, "command"):
        calc_dict.update({"command": calc.command})

    return calc_dict


def input2dict(atoms, calc=None, primitive=None, supercell=None, settings=None):
    """Convert metadata information to plain dict

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The structure to be converted to a dict
    calc: ase.calculators.calulator.Calculator
        The calculator to be converted to a dict
    primitive: ase.atoms.Atoms
        The primitive cell structure
    supercell: ase.atoms.Atoms
        The supercell cell structure
    settings: Settings
        The settings used to generate the inputs

    Returns
    -------
    input_dict: dict
        The dictionary representation of the inputs with the following items

        calculator: dict
            The calc_dict of the inputs
        atoms: dict
            The atoms_dict of the inputs
        primitive: dict
            The dict representation of the primitive cell structure
        supercell: dict
            The dict representation of the supercell cell structure
        settings: dict
            The dict representation of the settings used to generate the inputs
    """

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
        settings_dict = settings.to_dict()

        input_dict.update({"settings": settings_dict})

    return input_dict


def results2dict(atoms, calc=None):
    """extract information from atoms and calculator and convert to plain dict

    Args:
        atoms: The structure to be converted to a dict
        calc: The calculator to be converted to a dict

    Returns:
        dict: dictionary with items
            calculator: calc_dict of the inputs
            atoms: atoms_dict of the inputs
    """

    atoms_dict = atoms2dict(atoms)

    # calculated values
    calc_dict = {}

    if calc is None:
        calc = atoms.calc

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

    if not calc_dict:
        raise RuntimeError("calc_dict is empty, was the calculation successful?")

    return {"atoms": atoms_dict, "calculator": calc_dict}


def dict2atoms(atoms_dict, calc_dict=None, single_point_calc=True):
    """Convert dictionaries into atoms and calculator objects

    Parameters
    ----------
    atoms_dict: dict
        The dict representation of atoms
    calc_dict: dict
        The dict representation of calc

    Returns
    -------
    atoms: ase.atoms.Atoms
        atoms represented by atoms_dict with calculator represented by calc_dict
    """

    if "pbc" not in atoms_dict:
        atoms_dict.update({"pbc": "cell" in atoms_dict})

    uuid = atoms_dict.pop("unique_id", None)

    try:
        velocities = atoms_dict.pop("velocities")
    except KeyError:
        velocities = None

    if key_masses in atoms_dict:
        atoms_dict[key_masses] = expand_list(atoms_dict[key_masses])
    if key_symbols in atoms_dict:
        atoms_dict[key_symbols] = expand_list(atoms_dict[key_symbols])

    if key_constraints in atoms_dict:
        atoms_dict[key_constraints] = [
            dict2constraint(d) for d in atoms_dict[key_constraints]
        ]

    atoms = Atoms(**atoms_dict)

    if velocities is not None:
        atoms.set_velocities(velocities)

    # Calculator
    if calc_dict is not None:
        results = {}
        if "results" in calc_dict:
            results = calc_dict.pop("results")
        if single_point_calc:
            calc = SinglePointCalculator(atoms, **results)
            if "calculator" in calc_dict:
                calc.name = calc_dict["calculator"].lower()
            if "calculator_parameters" in calc_dict:
                calc.parameters.update(calc_dict["calculator_parameters"])
        else:
            calc = get_calculator_class(calc_dict["calculator"].lower())(
                **calc_dict["calculator_parameters"]
            )
            calc.results = results
        if "command" in calc_dict:
            calc.command = calc_dict["command"]
    else:
        calc = None

    atoms.calc = calc
    if "info" in atoms_dict:
        atoms.info = atoms_dict["info"]
    else:
        atoms.info = {}
    atoms.info["unique_id"] = uuid
    return atoms


def dict2json(dct, indent=0, outer=True):
    """convert python dictionary with scientific data to JSON

    Parameters
    ----------
    dct: dict
        Dictionary to convert to JSONAble format
    indent: int
        indentation for the json string
    outer: bool
        if True add outer { } to the string

    Returns
    -------
    rep: str
        json string of the dict
    """

    parts = []
    ind = indent * " "

    for key, val in dct.items():
        # convert np.ndarrays etc to lists
        if hasattr(val, "tolist") and callable(val.tolist):
            val = val.tolist()

        if isinstance(val, str):
            rep = f'"{val}"'
        elif isinstance(val, (float, np.float)):
            rep = "{1: .{0}e}".format(n_yaml_digits, val)
        elif isinstance(val, dict):
            # recursive formatting
            rep = f"{{\n{dict2json(val, 2*(1 + indent // 2), False)}}}"
        elif (
            isinstance(val, list)
            and len(list_dim(val)) == 2
            and list_dim(val)[1] == 3
            and isinstance(val[0][0], float)
        ):
            # this is most likely positions, velocities, forces, etc. -> format!
            rep = [
                " [{1: .{0}e}, {2: .{0}e}, {3: .{0}e}]".format(n_yaml_digits, *elem)
                for elem in val
            ]
            # join to have a comma separated list
            rep = f",\n{2*ind}".join(rep)
            # add leading [ and trailing ]
            rep = f"\n{2*ind}[{rep[1:]}"
            rep += "]"
        elif (
            isinstance(val, list)
            and len(list_dim(val)) == 3
            and list_dim(val)[1:3] == [3, 3]
            and isinstance(val[0][0][0], float)
        ):
            # this is most likely atomic stress -> format!
            rep = [
                "["
                + "[{1: .{0}e}, {2: .{0}e}, {3: .{0}e}]".format(n_yaml_digits, *elem[0])
                + f",\n{2*ind} "
                + "[{1: .{0}e}, {2: .{0}e}, {3: .{0}e}]".format(n_yaml_digits, *elem[1])
                + f",\n{2*ind} "
                + "[{1: .{0}e}, {2: .{0}e}, {3: .{0}e}]".format(n_yaml_digits, *elem[2])
                + "]"
                for elem in val
            ]
            # join to have a comma separated list
            rep = f",\n{2*ind}".join(rep)
            # add leading [ and trailing ]
            rep = f"\n{(2*ind)[:-1]}[{rep}"
            rep += "]"

        else:
            rep = json.dumps(val, cls=NumpyEncoder)

        parts.append(f'{ind}"{key}": {rep}')

    rep = ",\n".join(parts)

    if outer:
        rep = f"{{{rep}}}"
    # make sure only " are used to be understood by JSON
    return rep.replace("'", '"')


def get_json(obj):
    """Return json representation of obj

    Parameters
    ----------
    obj
        Object to convert to json

    Returns
    -------
    str
        json string of obj
    """
    return json.dumps(obj, cls=MyEncoder, sort_keys=True)


def atoms2json(atoms, reduce=True):
    """return json representation of Atoms"""
    rep = dict2json(atoms2dict(atoms, reduce=reduce))
    return rep


def json2atoms(rep):
    """return Atoms from json string"""
    atoms = dict2atoms(json.loads(rep))
    return atoms


def atoms_calc2json(
    atoms, ignore_results=False, ignore_keys=["unique_id"], ignore_calc_params=[]
):
    """Return json representation of atoms and calculator objects.

    Includes possibility to remove certain keys from the atoms dictionary,
    e.g. for hashing

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The structure to be converted to a json with attached calculator
    ignore_results: bool
        If True ignore the results in atoms.calc
    ignore_keys: list of str
        Ignore all keys in this list
    ignore_calc_params: list of str
        Ignore all keys in this list that represent calculator parameters

    Returns
    -------
    atoms_json: str
        Json representation of the atoms dictionary
    calc_json: str
        Json representation of the calculator dictionary
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