""" tools for storing MD trajectories

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""


# from vibes.helpers.converters import results2dict
# from vibes.helpers.converters import dict2json as dumper
from vibes.helpers.hash import hash_atoms

# from vibes import __version__ as version
# from vibes import son
from vibes.trajectory.io import metadata2file, reader, results2dict, step2file
from vibes.trajectory.trajectory import (
    Trajectory,
    key_reference_atoms,
    key_reference_positions,
    key_reference_primitive,
)

from .utils import get_hashes_from_trajectory_file
