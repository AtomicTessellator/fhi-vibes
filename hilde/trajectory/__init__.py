""" tools for storing MD trajectories 

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import numpy as np
from hilde.helpers.converters import input2dict, results2dict
from hilde.helpers.fileformats import to_yaml, from_yaml, last_from_yaml
from .reader import reader
