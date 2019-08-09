""" tools for storing MD trajectories

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

# from hilde import __version__ as version
# from hilde import son
# from hilde.helpers.converters import results2dict
# from hilde.helpers.converters import dict2json as dumper
from hilde.helpers import talk as _talk, Timer

_prefix = "trajectory"

Timer.prefix = _prefix


def talk(msg):
    """wrapper for `utils.talk` with prefix"""
    return _talk(msg, prefix=_prefix)


def get_hashes_from_trajectory(trajectory, verbose=False):
    """return all hashes from trajectory"""

    try:
        traj = reader(trajectory, verbose=verbose)
    except (FileNotFoundError, KeyError):
        return []

    hashes = []
    for atoms in traj:
        try:
            hashes.append(atoms.info["hash"])
        except (KeyError, AttributeError):
            hashes.append(hash_atoms(atoms))

    return hashes


from hilde.helpers.hash import hash_atoms
from hilde.trajectory.io import reader, metadata2file, step2file
from hilde.trajectory.trajectory import Trajectory, results2dict