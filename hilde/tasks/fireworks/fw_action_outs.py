from ase.symbols import Symbols
from ase.io.aims import read_aims
from fireworks import FWAction
from pathlib import Path
from phonopy import Phonopy
from phono3py.phonon3 import Phono3py
import shutil

import numpy as np

from hilde.helpers.converters import calc2dict, atoms2dict, dict2atoms
from hilde.helpers.fileformats import last_from_yaml
from hilde.helpers.k_grid import k2d
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import phonon_to_dict, phonon3_to_dict
from hilde.fireworks.workflow_generator import generate_firework

mod_name = __name__


def check_relaxation_complete(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that checks if a relaxation is converged (if outputs is True) and either stores the relaxed structure in the MongoDB or appends another Firework as its child to restart the relaxation
    Args:
        atoms: ASE Atoms object
            The original atoms at the start of this job
        calc: ASE Calculator object
            The original calculator
        outputs:
            The outputs from the function (assumes to be a single bool output)
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        func_kwargs: dict
            keyword arguments for func
        fw_settings: dict
            FireWorks specific settings
    """
    if "trajectory" in func_kwargs:
        last_step_dict = last_from_yaml(func_kwargs["trajectory"])
    elif "workdir" in func_kwargs:
        last_step_dict = last_from_yaml(
            Path(func_kwargs["workdir"] + "/trajectory.yaml").absolute()
        )
    else:
        last_step_dict = last_from_yaml(Path("./trajectory.yaml").absolute())

    for key, val in last_step_dict["atoms"].items():
        atoms[key] = val
    calc["results"] = last_step_dict["calculator"]
    for key, val in calc.items():
        atoms[key] = val
    next_step = last_step_dict["info"]["nsteps"] + 1

    if outputs:
        if "db_path" in func_fw_kwargs:
            db_path = func_fw_kwargs["db_path"]
            del (func_fw_kwargs["db_path"])
            update_phonon_db(db_path, atoms, atoms, **func_fw_kwargs)
        return FWAction(
            update_spec={
                fw_settings["out_spec_atoms"]: atoms,
                fw_settings["out_spec_calc"]: calc,
            }
        )
    del (calc["results"]["forces"])
    fw_settings["fw_name"] = fw_settings["fw_base_name"] + str(next_step)
    if "to_launchpad" in fw_settings and fw_settings["to_launchpad"]:
        fw_settings["to_launchpad"] = False
    if "trajectory" in func_kwargs:
        new_traj_list = func_kwargs["trajectory"].split(".")
    elif "workdir" in func_kwargs:
        new_traj_list = str(
            Path(func_kwargs["workdir"] + "/trajectory.yaml").absolute()
        ).split(".")
    else:
        new_traj_list = str(Path("." + "/trajectory.yaml").absolute()).split(".")

    try:
        temp_list = new_traj_list[-2].split("_")
        temp_list[-1] = str(int(temp_list[-1]) + 1)
        new_traj_list[-2] = "_".join(temp_list)
        func_kwargs["trajectory"] = ".".join(new_traj_list)
    except:
        new_traj_list[-2] += "_restart_1"
        func_kwargs["trajectory"] = ".".join(new_traj_list)

    fw = generate_firework(
        func=func,
        func_fw_out=func_fw_out,
        func_kwargs=func_kwargs,
        atoms=atoms,
        calc=calc,
        func_fw_out_kwargs=func_fw_kwargs,
        fw_settings=fw_settings,
    )
    return FWAction(detours=[fw])


def check_aims_relaxation_complete(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that checks if a relaxation is converged (if outputs is True) and either stores the relaxed structure in the MongoDB or appends another Firework as its child to restart the relaxation
    Args:
        atoms: ASE Atoms object
            The original atoms at the start of this job
        calc: ASE Calculator object
            The original calculator
        outputs:
            The outputs from the function (assumes to be a single bool output)
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        func_kwargs: dict
            keyword arguments for func
        fw_settings: dict
            FireWorks specific settings
    """
    func_fw_kwargs["relax_step"] += 1
    aims_out = open(func_kwargs["workdir"] + "/aims.out").readlines()
    converged = "Have a nice day" in aims_out[-2]
    calc = calc2dict(outputs.get_calculator())
    try:
        new_atoms = read_aims(func_kwargs["workdir"] + "/geometry.in.next_step")
        new_atoms.set_calculator(outputs.get_calculator())
        new_atoms = atoms2dict(new_atoms)
    except:
        if not converged:
            raise CalculationError(
                "There was a problem with the FHI Aims calculation stopping program here"
            )
        new_atoms = atoms2dict(outputs)
    for key, val in atoms['info'].items():
        print('\n', key, val, '\n')
        if key not in new_atoms['info']:
            new_atoms['info'][key] = val
    if converged:
        if "db_path" in func_fw_kwargs:
            db_path = func_fw_kwargs["db_path"]
            del (func_fw_kwargs["db_path"])
            update_phonon_db(db_path, new_atoms, new_atoms, **func_fw_kwargs)
        return FWAction(
            update_spec={
                fw_settings["out_spec_atoms"]: new_atoms,
                fw_settings["out_spec_calc"]: calc,
            }
        )

    del (calc["results"])
    fw_settings["fw_name"] = fw_settings["fw_base_name"] + str(
        func_fw_kwargs["relax_step"]
    )

    if "to_launchpad" in fw_settings and fw_settings["to_launchpad"]:
        fw_settings["to_launchpad"] = False

    fw = generate_firework(
        func=func,
        func_fw_out=func_fw_out,
        func_kwargs=func_kwargs,
        atoms=atoms,
        calc=calc,
        func_fw_out_kwargs=func_fw_kwargs,
        fw_settings=fw_settings,
    )
    return FWAction(detours=[fw])


def serial_phonopy_continue(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that checks if a set of force calculations are completed, and if not adds another Firework to the Workflow to continue the calculation
    Args:
        atoms: ASE Atoms object
            The original atoms at the start of this job
        calc: ASE Calculator object
            The original calculator
        outputs:
            The outputs from the function (assumes to be a single bool output)
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        func_kwargs: dict
            keyword arguments for func
        fw_settings: dict
            FireWorks specific settings
    """
    update_spec = {}
    if "in_spec_atoms" in fw_settings:
        update_spec[fw_settings["in_spec_atoms"]] = atoms
        update_spec[fw_settings["in_spec_calc"]] = calc
        at = fw_settings["in_spec_atoms"]
        cl = fw_settings["in_spec_calc"]
    else:
        at = atoms
        cl = calc
    if "kpoint_density_spec" in fw_settings:
        del(fw_settings["kpoint_density_spec"])

    if outputs:
        return FWAction(update_spec=update_spec)
    fw = generate_firework(
        func=func,
        func_fw_out=func_fw_out,
        func_kwargs=func_kwargs,
        atoms=at,
        calc=cl,
        func_fw_out_kwargs=func_fw_kwargs,
        fw_settings=fw_settings,
    )
    return FWAction(detours=[fw],update_spec=update_spec)


def check_kgrid_opt_completion(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that checks if an MD like calculation is converged (if outputs is True) and either stores the relaxed structure in the MongoDB or appends another Firework as its child to restart the MD
    Args:
        atoms: ASE Atoms object
            The original atoms at the start of this job
        calc: ASE Calculator object
            The original calculator
        outputs:
            The outputs from the function (assumes to be a single bool output)
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        func_kwargs: dict
            keyword arguments for func
        fw_settings: dict
            FireWorks specific settings
    """
    trajectory = Path(func_kwargs["workdir"]) / func_kwargs["trajectory"]
    last_step_dict = last_from_yaml(trajectory)
    for key, val in last_step_dict["atoms"].items():
        atoms[key] = val
    calc["results"] = last_step_dict["calculator"]
    for key, val in calc.items():
        atoms[key] = val
    next_step = last_step_dict["info"]["nsteps"] + 1
    if outputs[0]:
        up_spec = {
            fw_settings["out_spec_k_den"]: outputs[1],
            fw_settings["out_spec_atoms"]: atoms,
            fw_settings["out_spec_calc"]: calc2dict(outputs[2]),
        }
        return FWAction(update_spec=up_spec)

    fw_settings["fw_name"] = fw_settings["fw_base_name"] + str(next_step)
    if fw_settings["to_launchpad"]:
        fw_settings["to_launchpad"] = False
    new_traj_list = trajectory.split(".")
    try:
        temp_list = new_traj_list[-2].split("_")
        temp_list[-1] = str(int(temp_list[-1]) + 1)
        new_traj_list[-2] = "_".join(temp_list)
        trajectory = ".".join(new_traj_list)
    except:
        new_traj_list[-2] += "_restart_1"
        trajectory = ".".join(new_traj_list)

    fw = generate_firework(
        func=func,
        func_fw_out=func_fw_out,
        func_kwargs=func_kwargs,
        atoms=atoms,
        calc=outputs[2],
        func_fw_out_kwargs=func_fw_kwargs,
        fw_settings=fw_settings,
    )
    return FWAction(detours=[fw])


def fireworks_no_mods(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that does not change the FireWorks Workflow upon completion
    Args:
        atoms: ASE Atoms object
            The original atoms at the start of this job
        calc: ASE Calculator object
            The original calculator
        outputs:
            The outputs from the function (assumes to be a single bool output)
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        func_kwargs: dict
            keyword arguments for func
        fw_settings: dict
            FireWorks specific settings
    """
    return FWAction()


def fireworks_no_mods_gen_function(
    func, func_fw_out, *args, fw_settings=None, **kwargs
):
    """
    A function that does not change the FireWorks Workflow upon completion
    Args:
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        args: list
            List of arguments to pass to func
        fw_settings: dict
            FireWorks specific settings
        kwargs: dict
            keyword arguments for func
    """
    return FWAction()


def add_phonon_force_calcs(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that initializes and adds all phonopy force calculations to the FireWorks Workflow,
    and adds the Phonopy Object to the spec
    Args:
        atoms: ASE Atoms object
            The original atoms at the start of this job
        calc: ASE Calculator object
            The original calculator
        outputs:
            The outputs from the function (assumes to be a single bool output)
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        func_kwargs: dict
            keyword arguments for func
        fw_settings: dict
            FireWorks specific settings
    """
    detours = []
    update_spec = {}
    fw_settings = fw_settings.copy()
    if "spec" in fw_settings and "kpoint_density_spec" in fw_settings:
        fw_settings["spec"][fw_settings["kpoint_density_spec"]] = k2d(dict2atoms(atoms), calc["calculator_parameters"]["k_grid"])
    elif "kpoint_density_spec" in fw_settings:
        fw_settings["spec"] = {fw_settings["kpoint_density_spec"]: k2d(dict2atoms(atoms), calc["calculator_parameters"]["k_grid"])}
    if outputs[0]:
        update_spec["metadata_ph"] = outputs[0][4]
        calc_dict = calc2dict(outputs[0][0])
        fw_settings["mod_spec_add"] = "ph_forces"
        detours = add_to_detours(detours, func_fw_kwargs["phonopy_settings"], atoms, outputs[0][2], calc_dict, fw_settings, "ph")

    if outputs[1]:
        update_spec["metadata_ph3"] = outputs[1][4]
        calc_dict = calc2dict(outputs[1][0])
        fw_settings["mod_spec_add"] = "ph3_forces"
        detours = add_to_detours(detours, func_fw_kwargs["phono3py_settings"], atoms, outputs[1][2], calc_dict, fw_settings, "ph3")

    return FWAction(update_spec=update_spec, detours=detours)

def add_to_detours(detours, func_fw_kwargs, atoms, atoms_list, calc_dict, fw_settings, prefix):
    for i, sc in enumerate(atoms_list):
        if not sc:
            continue
        fw_settings=fw_settings.copy()
        sc.info["displacement_id"] = i
        sc_dict = atoms2dict(sc)
        for key, val in calc_dict.items():
            sc_dict[key] = val
        calc_kwargs = {"workdir": func_fw_kwargs["workdir"] + f"/{i:05d}"}
        fw_settings["fw_name"] = prefix + f"forces_{Symbols(atoms['numbers']).get_chemical_formula()}_{i}"
        detours.append(
            generate_firework(
                func="hilde.tasks.calculate.calculate",
                func_fw_out="hilde.tasks.fireworks.fw_action_outs.mod_spec_add",
                func_kwargs=calc_kwargs,
                atoms=sc_dict,
                calc=calc_dict,
                atoms_calc_from_spec=False,
                fw_settings=fw_settings,
            )
        )
    return detours

def mod_spec_add(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that appends the current results to a specified spec in the MongoDB
    Args:
        atoms: ASE Atoms object
            The original atoms at the start of this job
        calc: ASE Calculator object
            The original calculator
        outputs:
            The outputs from the function (assumes to be a single bool output)
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        func_kwargs: dict
            keyword arguments for func
        fw_settings: dict
            FireWorks specific settings
    """
    atoms_dict = atoms2dict(outputs)
    return FWAction(mod_spec=[{"_push": {fw_settings["mod_spec_add"]: atoms_dict}}])
