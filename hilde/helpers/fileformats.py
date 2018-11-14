""" tools for conerting atoms objects to json representations """


import json
from pathlib import Path
import yaml
import numpy as np

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


def list2str(lis):
    """convert list to string"""
    return "[{}]".format(", ".join([str(el) for el in lis]))


class NumpyEncoder(json.JSONEncoder):
    """ Decode numerical objects that json cannot parse by default"""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        if isinstance(obj, complex):
            return (float(obj.real), float(obj.imag))
        return super().default(self, obj)


def to_yaml(obj, file, mode="a"):
    """ Dump a python object ot file """

    # backup
    if mode == "w":
        if Path(file).exists():
            Path(file).rename(f"{file}.bak")

    with open(file, mode) as f:
        yaml.dump(obj, f)


def from_yaml(file):
    """ return from yaml file """
    with Path(file).open() as f:
        return yaml.load(f, Loader=Loader)


def from_json(file):
    """ return from json file """

    with Path(file).open() as f:
        return json.load(f)


def to_json(obj, f, indent=1):
    """ write array (or similar) to json file """

    with Path(f).open("w") as f:
        json.dump(obj, f, cls=NumpyEncoder, indent=indent)


def append_to_json_array(obj, f, indent=1):
    """ append contents to existing json array """

    if not Path(f).exists():
        to_json([obj], f, indent=indent)
        return

    new_array = from_json(f)
    new_array.append(obj)

    to_json(new_array, f, indent=indent)
