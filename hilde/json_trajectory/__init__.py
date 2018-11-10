""" A trajectory for e.g. MD simulations using JSON """

import json
from pathlib import Path
import numpy as np


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


def array_from_json(f):
    """ return array from json file """

    with Path(f).open() as f:
        return json.load(f)


def array_to_json(array, f, indent=2):
    """ write array to json file """

    with Path(f).open("w") as f:
        json.dump(array, f, cls=NumpyEncoder, indent=indent)


def append_to_json_array(obj, f):
    """ append contents of array to existing json file """

    if not Path(f).exists():
        array_to_json([obj], f)
        return

    new_array = array_from_json(f)

    new_array.append(obj)

    array_to_json(new_array, f)
