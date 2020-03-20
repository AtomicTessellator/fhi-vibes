"""Post processing for FHI-aims calculations"""

from pathlib import Path
from shutil import copyfile

import numpy as np
from ase.io.aims import read_aims

from vibes.filenames import filenames
from vibes.fireworks.tasks.calculate_wrapper import check_if_failure_ok
from vibes.helpers.converters import atoms2dict, calc2dict, dict2atoms, key_constraints


def check_aims(atoms, calculator, outputs, **kwargs):
    """
    A function that checks if a relaxation is converged (if outputs is True) and either
    stores the relaxed structure in the MongoDB or appends another Firework as its child
    to restart the relaxation
    Args:
        atoms (ASE Atoms object): The original atoms at the start of this job
        calculator (ASE Calculator object): The original calculator
        outputs (ASE Atoms Object): The geometry of the final relaxation step
    Returns (FWAction): The correct action if convergence is reached
    """
    calc_number = kwargs.get("calc_number", 0) + 1
    path = kwargs["workdir"]
    aims_out = np.array(open(path / filenames.output.aims).readlines())
    completed = "Have a nice day" in aims_out[-2] or "Have a nice day" in aims_out[-3]
    calculator = calc2dict(outputs.get_calculator())
    walltime = kwargs.get("walltime", 0)
    try:
        if "relax_geometry" in calculator["calculator_parameters"]:
            workdir = Path(kwargs["workdir"])
            new_atoms = read_aims(workdir / filenames.atoms_next)
            new_atoms.set_calculator(outputs.get_calculator())
            new_atoms.info = atoms["info"].copy()
        else:
            new_atoms = dict2atoms(atoms)
    except FileNotFoundError:
        if not completed:
            failure_ok = check_if_failure_ok(aims_out, walltime)
            if failure_ok:
                walltime *= 2
                calculator.parameters["walltime"] = walltime
            else:
                raise IOError(
                    "There was a problem with the FHI Aims calculation stopping here"
                )
        new_atoms = outputs
    new_atoms_dict = atoms2dict(new_atoms)
    new_atoms_dict[key_constraints] = atoms.get(key_constraints, ())
    for key, val in atoms["info"].items():
        if key not in new_atoms_dict["info"]:
            new_atoms_dict["info"][key] = val
    path = Path(kwargs["workdir"]) / filenames.output.aims
    copyfile(f"{path}", f"{path}.{calc_number}")
    return completed, calc_number, new_atoms_dict, walltime
