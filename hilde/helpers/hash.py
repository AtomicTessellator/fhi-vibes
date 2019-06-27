""" Tools for hashing atoms objects """

from json import dumps
from pathlib import Path
from configparser import ConfigParser
from hashlib import sha1 as hash_sha
from .converters import atoms2dict, atoms2json, get_json, dict2json


def hashfunc(string, empty_str="", digest=True):
    """Wrap the sha hash function and check for empty objects

    Parameters
    ----------
    string: str
        string to hash
    emptystr: str
        What an empy string should map to
    digest: bool
        If True digest the hash

    Returns
    -------
    str
        Hash of string
    """
    if string in ("", "[]", "{}", "None"):
        string = empty_str
    if digest:
        return hash_sha(string.encode("utf8")).hexdigest()
    return hash_sha(string.encode("utf8"))


def hash_atoms(atoms):
    """hash the atoms object as it would be written to trajectory

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Atoms to has

    Returns
    -------
    atoms_hash: str
        The hash of the atoms Object
    """
    a = atoms.copy()
    a.info = {}

    rep = dict2json(atoms2dict(a))

    atoms_hash = hashfunc(rep)

    return atoms_hash


def hash_atoms_and_calc(
    atoms,
    ignore_results=True,
    ignore_keys=["unique_id", "info"],
    ignore_calc_params=[],
    ignore_file=None,
):
    """Hash atoms and calculator object, with possible ignores

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
    ignore_file: str
        Path to a file with standard keys to ignore

    Returns
    -------
    atomshash: str
        hash of the atoms
    calchash: str
        hash of atoms.calc
    """

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

    atomshash = hashfunc(atomsjson)
    calchash = hashfunc(calcjson)

    return atomshash, calchash


def hash_traj(ca, meta, hash_meta=False):
    """hash of a trajectory file

    Parameters
    ----------
    ca: list of ase.atoms.Atoms
        Atoms objects inside a trajectory file
    meta: dict
        Metadata for the trajectory
    hash_meta: bool
        if True hash the meta data

    Returns
    atomshash: str
        hash of all of the ca objects
    metahash: str
        hash of the metadata
    """

    ca_dct = [atoms2json(at) for at in ca]
    dct = dict(meta, calculated_atoms=ca_dct)
    if hash_meta:
        return hashfunc(dumps(dct)), hashfunc(dumps(meta))
    return hashfunc(dumps(dct))


def hash_dict(dct):
    """hash a dictionary and check if for species_dir is a key, if so remove it

    Parameters
    ----------
    dct: dict
        Dictionary to hash

    Returns
    -------
    str
        hash of the dictionary
    """
    if "calculator_parameters" in dct:
        if "species_dir" in dct["calculator_parameters"]:
            dct["calculator_parameters"]["species_dir"] = Path(
                dct["calculator_parameters"]["species_dir"]
            ).parts[-1]

    return hashfunc(get_json(dct))
