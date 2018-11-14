""" tools for storing MD trajectories """

import numpy as np
from ase import units as u
from hilde.helpers.fileformats import to_yaml, from_yaml
from .reader import reader


def step2file(atoms, calc, md, file="md_trajectory.yaml"):
    """ Save the current state of MD to file """

    to_yaml([step2dict(atoms, calc, md)], file)


def metadata2file(atoms, calc, md, file="md_metadata.yaml"):
    """ save MD metadata to file """

    metadata = metadata2dict(atoms, calc, md)

    to_yaml([metadata], file, mode="w")


def step2dict(atoms, calc, md):
    """ extract information from md step and convet to plain dict """

    # info from MD algorithm
    md_dict = {"nsteps": md.nsteps, "dt": md.dt}

    # structure
    atoms_dict = {
        "numbers": atoms.numbers.tolist(),
        "positions": atoms.positions.tolist(),
    }

    if atoms.get_velocities() is not None:
        atoms_dict.update({"velocities": atoms.get_velocities().tolist()})

    # if periodic system, append lattice
    if any(atoms.pbc):
        atoms_dict.update({"cell": atoms.cell.tolist()})

    # calculated values
    calc_dict = {}
    # convert numpy arrays into ordinary lists
    for key, val in calc.results.items():
        if isinstance(val, np.ndarray):
            calc_dict[key] = val.tolist()
        elif isinstance(val, np.float):
            calc_dict[key] = float(val)
        else:
            calc_dict[key] = val

    return {"MD": md_dict, "atoms": atoms_dict, "calculator": calc_dict}


def metadata2dict(atoms, calc, md):
    """ convert metadata information to plain dict """

    md_dict = md.todict()
    # save time and mass unit
    md_dict.update({"fs": u.fs, "dt": md.dt, "kg": u.kg})

    # structure
    atoms_dict = {
        "numbers": atoms.numbers.tolist(),
        "symbols": [f"{sym}" for sym in atoms.symbols],
        "masses": md.masses.T.tolist()[0],
        "positions": atoms.positions.tolist(),
    }
    # if periodic system, append lattice
    if any(atoms.pbc):
        atoms_dict.update({"cell": atoms.cell.tolist()})

    params = calc.todict()
    for key, val in params.items():
        if isinstance(val, tuple):
            params['key'] = list(val)

    calc_dict = {"name": calc.__class__.__name__, "params": params}

    return {"MD": md_dict, "atoms": atoms_dict, "calculator": calc_dict}
