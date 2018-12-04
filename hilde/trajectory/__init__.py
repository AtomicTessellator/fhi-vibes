""" tools for storing MD trajectories 

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import numpy as np
from hilde.helpers.converters import input2dict, results2dict
from hilde.helpers.fileformats import to_yaml, from_yaml, last_from_yaml
from hilde.helpers.hash import hash_atoms
from .reader import reader


def get_hashes_from_trajectory(trajectory):
    """ return all hashes from trajectory """

    try:
        traj = reader(trajectory)
    except (FileNotFoundError, KeyError):
        return []

    hashes = []
    for atoms in traj:
        try:
            hashes.append(atoms.info["hash"])
        except AttributeError:
            hashes.append(hash_atoms(atoms))

    return hashes
