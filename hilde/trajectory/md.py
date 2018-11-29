""" tools for storing MD trajectories 

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import numpy as np
from ase import units as u
from . import input2dict, results2dict
from hilde.helpers.fileformats import to_yaml, from_yaml, last_from_yaml


def step2file(atoms, calc, md, file="md_trajectory.yaml", NPT=False):
    """ Save the current state of MD to file """

    to_yaml(step2dict(atoms, calc, md, NPT), file)


def metadata2file(atoms, calc, md, file="md_metadata.yaml"):
    """ save MD metadata to file """

    metadata = metadata2dict(atoms, calc, md)

    to_yaml(metadata, file, mode="w")


def step2dict(atoms, calc, md, NPT=False):
    """ extract information from md step and convet to plain dict """

    # info from MD algorithm
    md_dict = {"nsteps": md.nsteps, "dt": md.dt}

    return {"MD": md_dict, **results2dict(atoms, calc, append_cell=NPT)}


def metadata2dict(atoms, calc, md):
    """ convert metadata information to plain dict """

    md_dict = md.todict()
    # save time and mass unit
    md_dict.update({"fs": u.fs, "kB": u.kB, "dt": md.dt, "kg": u.kg})

    return {"MD": md_dict, **input2dict(atoms, calc)}
