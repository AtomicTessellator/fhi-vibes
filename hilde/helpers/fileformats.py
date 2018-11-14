""" tools for conerting atoms objects to json representations """


import json
from pathlib import Path
import yaml
import numpy as np

try:
    from yaml import CSafeLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import SafeLoader as Loader, Dumper


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
            if "traj" in file:
                raise Exception("Possibly overwriting a trajectory, please check")
            Path(file).rename(f"{file}.bak")

    with open(file, mode) as f:
        if "a" in mode:
            f.write("---\n")
        yaml.dump(obj, f)


def from_yaml(file):
    """ return from yaml file """
    with Path(file).open() as f:
        stream = yaml.load_all(f, Loader=Loader)
        return [elem for elem in stream]


def from_json(file):
    """ return from json file """

    blobs = Path(file).read_text().split('---')

    return [json.loads(blob) for blob in blobs if blob.strip()]


def to_json(obj, file, mode="a", indent=1):
    """ Dump a python object to json file """

    # backup
    if mode == "w":
        if Path(file).exists():
            if "traj" in file:
                raise Exception("Possibly overwriting a trajectory, please check")
            Path(file).rename(f"{file}.bak")

    with open(file, mode) as f:
        if "a" in mode:
            f.write("\n---\n")
        json.dump(obj, f, cls=NumpyEncoder, indent=indent)


def append_to_json_array(obj, file, indent=1):
    """ append contents to existing json array """

    if not Path(file).exists():
        to_json([obj], file, indent=indent)
        return

    new_array = from_json(file)
    new_array.append(obj)

    to_json(new_array, file, indent=indent)
