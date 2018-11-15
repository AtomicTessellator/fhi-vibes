""" tools for storing MD trajectories 

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

from hilde.helpers.fileformats import to_yaml, from_yaml, last_from_yaml
from .reader import reader
from .md import step2file, metadata2file
