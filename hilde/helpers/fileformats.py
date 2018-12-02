""" tools for conerting atoms objects to json representations """


import json
from pathlib import Path
import yaml
import numpy as np

try:
    from yaml import CSafeLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import SafeLoader as Loader, Dumper

from hilde.helpers import list_dim
from hilde.konstanten.io import n_yaml_digits


class NumpyEncoder(json.JSONEncoder):
    """ Decode numerical objects that json cannot parse by default"""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        if isinstance(obj, complex):
            return (float(obj.real), float(obj.imag))
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)


def backup(file):
    if Path(file).exists():
        if "traj" in str(file):
            raise Exception("Possibly overwriting a trajectory, please check")
        Path(file).rename(f"{file}.bak")


def to_yaml(obj, file, mode="a", use_json=True):
    """ Dump a python object ot file """
    if use_json:
        to_json(obj, file, mode)
        return

    # backup
    if mode == "w":
        backup(file)

    with open(file, mode) as f:
        if "a" in mode:
            f.write("---\n")
        yaml.dump(obj, f)


def from_yaml(file, use_json=True):
    """ return from yaml file """
    with Path(file).open() as f:
        if use_json:
            return from_json(file)
        stream = yaml.load_all(f, Loader=Loader)
        return [elem for elem in stream]


def last_from_yaml(file):
    """ return last entry from yaml file """

    blobs = Path(file).read_text().split("---")

    return yaml.load(blobs[-1], Loader=Loader)


def dict2json(dct, indent=0, outer=True):
    """ convert pyton dictionary with scientific data to JSON """

    parts = []
    ind = indent * " "

    for key, val in dct.items():
        if isinstance(val, str):
            rep = f'"{val}"'
        elif isinstance(val, dict):
            # recursive formatting
            rep = f"{{\n{dict2json(val, 2*(1 + indent // 2), False)}}}"
        elif (
            isinstance(val, list)
            and len(list_dim(val)) == 2
            and type(val[0][0]) == float
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
        else:
            rep = json.dumps(val, cls=NumpyEncoder)

        parts.append(f'{ind}"{key}": {rep}')

    rep = ",\n".join(parts)

    if outer:
        rep = f"{{{rep}}}"
    # make sure only " are used to be understood by JSON
    return rep.replace("'", '"')


def from_json(file):
    """ return from json file """

    blobs = Path(file).read_text().split("---")

    return [json.loads(blob) for blob in blobs if blob.strip()]


def to_json(obj, file, mode="a", indent=1):
    """ Dump a python object to json file """

    # backup
    if mode == "w":
        backup(file)

    if isinstance(obj, dict):
        rep = dict2json(obj)
    elif isinstance(obj, list) and isinstance(obj[0], dict):
        reps = [dict2json(elem) for elem in obj]
        rep = "[" + ",\n".join(reps) + "]"
    else:
        rep = json.dumps(obj, cls=NumpyEncoder, indent=indent)

    with open(file, mode) as f:
        if "a" in mode:
            f.write("\n---\n")
        f.write(rep)
