""" tools for storing MD trajectories 

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import numpy as np
from hilde.helpers.fileformats import to_yaml, from_yaml, last_from_yaml
from .reader import reader


def input2dict(atoms, calc=None):
    """ convert metadata information to plain dict """

    if calc is None:
        calc = atoms.calc

    # structure
    atoms_dict = {
        "symbols": [f"{sym}" for sym in atoms.symbols],
        "masses": atoms.get_masses().tolist(),
        "positions": atoms.positions.tolist(),
    }

    # if periodic system, append lattice
    if any(atoms.pbc):
        atoms_dict.update({"cell": atoms.cell.tolist()})

    if calc is None:
        return {"atoms": atoms_dict, "calculator": {}}

    params = calc.todict()
    for key, val in params.items():
        if isinstance(val, tuple):
            params[key] = list(val)

    calc_dict = {"name": calc.__class__.__name__, "params": params}

    return {"atoms": atoms_dict, "calculator": calc_dict}


def results2dict(atoms, calc, append_cell=False):
    """ extract information from atoms and calculator and convert to plain dict """

    # structure
    atoms_dict = {"positions": atoms.positions.tolist()}

    if atoms.get_velocities() is not None:
        atoms_dict.update({"velocities": atoms.get_velocities().tolist()})

    # if periodic system, append lattice
    if append_cell and any(atoms.pbc):
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

    return {"atoms": atoms_dict, "calculator": calc_dict}
