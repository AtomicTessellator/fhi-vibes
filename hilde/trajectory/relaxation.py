""" tools for storing MD trajectories 

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import numpy as np
from ase import units as u
from . import input2dict, results2dict
from hilde.helpers.fileformats import to_yaml, from_yaml, last_from_yaml


def step2file(atoms, calc, opt, file="relax_trajectory.yaml", unit_cell=True):
    """ Save the current state of MD to file """

    to_yaml(step2dict(atoms, calc, opt, unit_cell=unit_cell), file)


def metadata2file(atoms, calc, opt, file="opt_metadata.yaml"):
    """ save opt metadata to file """

    metadata = metadata2dict(atoms, calc, opt)

    to_yaml(metadata, file, mode="w")


def step2dict(atoms, calc, opt, unit_cell=True):
    """ extract information from opt step and convet to plain dict """

    # info from opt algorithm
    opt_dict = {"nsteps": opt.nsteps}

    return {"opt": opt_dict, **results2dict(atoms, calc, append_cell=unit_cell)}


def metadata2dict(atoms, calc, opt):
    """ convert metadata information to plain dict """

    opt_dict = opt.todict()

    return {"opt": opt_dict, **input2dict(atoms, calc)}
