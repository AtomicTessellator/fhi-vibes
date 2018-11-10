""" Tools for hashing atoms objects """

from pathlib import Path
from configparser import ConfigParser
from hashlib import sha1 as hash_sha
from .fileformats import atoms2json


def hashfunc(string, empty_str=""):
    """ Wrap the sha hash function and check for empty objects """
    if string in ("", "[]", "{}", "None"):
        string = empty_str
    return hash_sha(string.encode("utf8"))


def hash_atoms(
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

            ignore_keys += [
                key
                for key in configparser["atoms"]
                if not configparser.getboolean("atoms", key)
            ]

            ignore_calc_params = [
                key
                for key in configparser["calculator"]
                if not configparser.getboolean("calculator", key)
            ]

    atomsjson, calcjson = atoms2json(
        atoms, ignore_results, ignore_keys, ignore_calc_params
    )

    atomshash = hashfunc(atomsjson).hexdigest()
    calchash = hashfunc(calcjson).hexdigest()

    return atomshash, calchash
