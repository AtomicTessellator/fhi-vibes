"""Wrap son to provide pretty formatted json"""

import son

from vibes.helpers import warn
from vibes.helpers.converters import dict2json


def dump(*args, **kwargs):
    """Wrapper for son.dump"""
    return son.dump(*args, **{"dumper": dict2json, **kwargs})


def load(*args, **kwargs):
    """Wrapper for son.load"""
    return son.load(*args, **kwargs)


def open(*args, **kwargs):
    """Wrapper for son.open"""
    return son.open(*args, **kwargs)


def last_from(file, allow_empty=False):
    """
    Return last entry from son file

    Parameters
    ----------
    file: str
        Path to file to load

    Returns
    -------
    data[-1]: dict
        Last entry in the son file

    """
    _, data = son.load_last(file)
    if not allow_empty and data is None:
        warn(
            "** trajectory lacking the first step, please CHECK!",
            level=2,
        )

    return data
