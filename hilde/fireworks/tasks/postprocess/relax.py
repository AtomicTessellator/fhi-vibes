"""Post processing for FHI-aims calculations"""

from shutil import copyfile

from ase.io.aims import read_aims

import numpy as np

from hilde.fireworks.tasks.calculate_wrapper import check_if_failure_ok
from hilde.helpers.converters import atoms2dict, dict2atoms, calc2dict, key_constraints


def check_aims(atoms, calc, outputs, **kwargs):
    """
    A function that checks if a relaxation is converged (if outputs is True) and either
    stores the relaxed structure in the MongoDB or appends another Firework as its child
    to restart the relaxation
    Args:
        atoms (ASE Atoms object): The original atoms at the start of this job
        calc (ASE Calculator object): The original calculator
        outputs (ASE Atoms Object): The geometry of the final relaxation step
    Returns (FWAction): The correct action (restart or updated spec) if convergence is reached
    """
    calc_number = kwargs.get("calc_number", 0) + 1
    aims_out = np.array(open(kwargs["workdir"] + "/aims.out").readlines())
    completed = "Have a nice day" in aims_out[-2] or "Have a nice day" in aims_out[-3]
    calc = calc2dict(outputs.get_calculator())
    walltime = kwargs.get("walltime", 0)
    try:
        if "relax_geometry" in calc["calculator_parameters"]:
            new_atoms = read_aims(kwargs["workdir"] + "/geometry.in.next_step")
            new_atoms.set_calculator(outputs.get_calculator())
            new_atoms.info = atoms["info"].copy()
        else:
            new_atoms = dict2atoms(atoms)
    except FileNotFoundError:
        if not completed:
            failure_ok = check_if_failure_ok(aims_out, walltime)
            if failure_ok:
                walltime *= 2
                calc.parameters["walltime"] = walltime
            else:
                raise IOError(
                    "There was a problem with the FHI Aims calculation stopping program here"
                )
        new_atoms = outputs
    new_atoms_dict = atoms2dict(new_atoms)
    new_atoms_dict[key_constraints] = atoms.get(key_constraints, list())
    for key, val in atoms["info"].items():
        if key not in new_atoms_dict["info"]:
            new_atoms_dict["info"][key] = val
    copyfile(
        f"{kwargs['workdir']}/aims.out", f"{kwargs['workdir']}/aims.out.{calc_number}"
    )
    return completed, calc_number, new_atoms_dict, walltime
