import os
from pathlib import Path

import numpy as np

from fireworks import FWAction
from ase.io.aims import read_aims

from hilde.fireworks.workflows.firework_generator import generate_firework, time2str
from hilde.phonon_db.ase_converters import atoms2dict, calc2dict
from hilde.helpers.fileformats import last_from_yaml
from hilde.helpers.k_grid import k2d
from hilde.helpers.watchdogs import str2time

def check_aims(
    atoms, calc, outputs, **kwargs
):
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
    completed = "Have a nice day" in aims_out[-2]
    calc = calc2dict(outputs.get_calculator())
    walltime = kwargs.get("walltime", 0)
    try:
        if "relax_geometry" in calc["calculator_parameters"]:
            new_atoms = read_aims(kwargs["workdir"] + "/geometry.in.next_step")
            new_atoms.set_calculator(outputs.get_calculator())
            new_atoms.info = atoms["info"].copy()
        else:
            new_atoms = atoms.copy()
    except FileNotFoundError:
        if not completed:
            line_sum = np.where(
                aims_out
                == "          Detailed time accounting                     :  max(cpu_time)    wall_clock(cpu1)\n"
            )[0]
            sum_present = len(line_sum) > 0
            time_used = float(aims_out[line_sum[0] + 1].split(":")[1].split("s")[1])

            if (sum_present and time_used / walltime > 0.95):
                walltime *= 2
                calc.parameters["walltime"] = walltime
            else:
                raise IOError(
                    "There was a problem with the FHI Aims calculation stopping program here"
                )
        new_atoms = outputs
    new_atoms_dict = atoms2dict(new_atoms)
    for key, val in atoms["info"].items():
        if key not in new_atoms_dict["info"]:
            new_atoms_dict["info"][key] = val
    Path(f"{kwargs['workdir']}/aims.out").rename(f"{kwargs['workdir']}/aims.out.{calc_number}")
    return completed, calc_number, new_atoms_dict, walltime