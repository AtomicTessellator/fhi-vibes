"""Wrap son to provide pretty formatted json"""

import son
from hilde.helpers.converters import dict2json


def dump(*args, **kwargs):
    """wrapper for son.dump"""
    return son.dump(*args, **{"dumper": dict2json, **kwargs})


def load(*args, **kwargs):
    """wrapper for son.load"""
    return son.load(*args, **kwargs)

def last_from(file):
    """ return last entry from yaml file """

    _, data = son.load(file)

    return data[-1]
