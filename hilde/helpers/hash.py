from hashlib import sha1 as hash_sha
from ase.db.row import atoms2dict
from ase.io.jsonio import encode, MyEncoder
from ase.calculators.calculator import all_properties
from pathlib import Path
from configparser import ConfigParser
import json


def get_json(obj):
    "Return json representation of obj"
    return json.dumps(obj, cls=MyEncoder, sort_keys=True)


def hashfunc(string, default=""):
    """ Wrap the sha hash function and check for empty objects """
    if string == "" or string == "[]" or string == "{}" or string == "None":
        string = ""
    return hash_sha(string.encode("utf8"))


def atoms2json(atoms, ignore_atoms=[], ignore_calc=[]):
    atomsdict = atoms2dict(atoms)

    for name in ["unique_id"] + ignore_atoms:
        atomsdict.pop(name)
    calcdict = {}

    # clean calculator entries
    if "calculator_parameters" in atomsdict:
        calculator_params = atomsdict["calculator_parameters"]
        for key in calculator_params:
            if key in ignore_calc:
                calculator_params.pop(key)

        if "species_dir" in calculator_params:
            calculator_params["species_dir"] = Path(
                calculator_params["species_dir"]
            ).parts[-1]

    for name in all_properties:
        atomsdict.pop(name, None)

    for name in ["calculator", "calculator_parameters"]:
        if name in atomsdict:
            calcdict[name] = atomsdict.pop(name)

    return get_json(atomsdict), get_json(calcdict)


def hash_atoms(atoms, ignore_atoms=[], ignore_calc=[], ignore_file=None):
    """ Hash atoms and calculator object, with possible ignores"""
    if ignore_file is not None:
        fil = Path(ignore_file)
        if fil.exists():
            configparser = ConfigParser()
            configparser.read(fil)

            ignore_atoms = [
                key
                for key in configparser["atoms"]
                if not configparser.getboolean("atoms", key)
            ]

            ignore_calc = [
                key
                for key in configparser["calculator"]
                if not configparser.getboolean("calculator", key)
            ]

    atomsjson, calcjson = atoms2json(atoms, ignore_atoms, ignore_calc)

    atomshash = hashfunc(atomsjson).hexdigest()
    calchash = hashfunc(calcjson).hexdigest()

    return atomshash, calchash
