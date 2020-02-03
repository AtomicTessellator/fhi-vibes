""" tools for storing MD trajectories

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

# from vibes import __version__ as version
# from vibes import son
# from vibes.helpers.converters import results2dict
# from vibes.helpers.converters import dict2json as dumper
from vibes.helpers import Timer
from vibes.helpers import talk as _talk
from vibes.helpers.hash import hash_atoms
from vibes.trajectory.io import metadata2file, reader, results2dict, step2file
from vibes.trajectory.trajectory import Trajectory

_prefix = "trajectory"
_fc_key = "force_constants"

key_reference_atoms = "reference atoms"
key_reference_primitive = "reference primitive atoms"
key_reference_positions = "reference positions"


Timer.prefix = _prefix


def talk(msg):
    """wrapper for `utils.talk` with prefix"""
    return _talk(msg, prefix=_prefix)


def get_hashes_from_trajectory_file(trajectory, verbose=False):
    """return all hashes from trajectory"""

    try:
        traj = reader(trajectory, verbose=verbose)
    except (FileNotFoundError, KeyError):
        return []

    return traj.get_hashes()
