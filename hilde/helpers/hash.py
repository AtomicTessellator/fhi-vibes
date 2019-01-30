""" Tools for hashing atoms objects """

from json import dumps
from pathlib import Path
from configparser import ConfigParser
from hashlib import sha1 as hash_sha
from .converters import atoms2json

def hashfunc(string, empty_str=""):
    """ Wrap the sha hash function and check for empty objects """
    if string in ("", "[]", "{}", "None"):
        string = empty_str
    return hash_sha(string.encode("utf8"))


def hash_atoms_and_calc(
    atoms,
    ignore_results=True,
    ignore_keys=["unique_id", "info"],
    ignore_calc_params=[],
    ignore_file=None,
):
    """ Hash atoms and calculator object, with possible ignores"""

    if ignore_file is not None:
        fil = Path(ignore_file)
        if fil.exists():
            configparser = ConfigParser()
            configparser.read(fil)
            ignores = configparser["hash_ignore"]

            ignore_keys += [key for key in ignores if not ignores.getboolean(key)]

            ignore_calc_params = [key for key in ignores if not ignores.getboolean(key)]

    atomsjson, calcjson = atoms2json(
        atoms, ignore_results, ignore_keys, ignore_calc_params
    )

    atomshash = hashfunc(atomsjson).hexdigest()
    calchash = hashfunc(calcjson).hexdigest()

    return atomshash, calchash

def hash_traj(ca, meta, hash_meta=False):
    ca_dct = [atoms2json(at) for at in ca]
    dct = dict(meta, calculated_atoms=ca_dct)
    if hash_meta:
        return hashfunc(dumps(dct)).hexdigest(), hashfunc(dumps(meta)).hexdigest()
    return hashfunc(dumps(dct)).hexdigest()

def hash_atoms(atoms):
    """ hash only the atoms object """
    hash_atoms = atoms.copy()
    hash_atoms.info = {}

    atoms_json, _ = atoms2json(atoms)

    atoms_hash = hashfunc(atoms_json).hexdigest()

    return atoms_hash
