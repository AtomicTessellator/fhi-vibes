""" A trajectory for e.g. MD simulations using JSON """

import re
import json
from pathlib import Path
import numpy as np


def list2str(lis):
    return "[{}]".format(", ".join([str(el) for el in lis]))


class NoIndent(object):
    """ Value wrapper. """

    def __init__(self, value):
        self.value = value


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


class MyEncoder(NumpyEncoder):
    FORMAT_SPEC = "@@{}@@"
    regex = re.compile(FORMAT_SPEC.format(r"(\d+)"))

    def __init__(self, **kwargs):
        # Save copy of any keyword argument values needed for use here.
        self.__sort_keys = kwargs.get("sort_keys", None)
        super(MyEncoder, self).__init__(**kwargs)

    def default(self, obj):
        return (
            self.FORMAT_SPEC.format(id(obj))
            if isinstance(obj, NoIndent)
            else super(MyEncoder, self).default(obj)
        )

    def encode(self, obj):
        format_spec = self.FORMAT_SPEC  # Local var to expedite access.
        json_repr = super(MyEncoder, self).encode(obj)  # Default JSON.

        # Replace any marked-up object ids in the JSON repr with the
        # value returned from the json.dumps() of the corresponding
        # wrapped Python object.
        for match in self.regex.finditer(json_repr):
            # see https://stackoverflow.com/a/15012814/355230
            id = int(match.group(1))
            no_indent = PyObj_FromPtr(id)
            json_obj_repr = json.dumps(no_indent.value, sort_keys=self.__sort_keys)

            # Replace the matched id string with json formatted representation
            # of the corresponding Python object.
            json_repr = json_repr.replace(
                '"{}"'.format(format_spec.format(id)), json_obj_repr
            )

        return json_repr


def from_json(f):
    """ return from json file """

    with Path(f).open() as f:
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
