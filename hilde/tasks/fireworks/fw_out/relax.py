"""FWAction generators for relaxations"""
from pathlib import Path

from fireworks import FWAction
from ase.io.aims import read_aims

from hilde.fireworks.workflows.workflow_generator import generate_firework
from hilde.helpers.converters import atoms2dict, calc2dict
from hilde.helpers.fileformats import last_from_yaml
from hilde.helpers.k_grid import k2d


def check_relaxation_complete(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that checks if a relaxation is converged (if outputs is True) and either
    stores the relaxed structure in the MongoDB or appends another Firework as its child
    to restart the relaxation
    Args:
        atoms (ASE Atoms object): The original atoms at the start of this job
        calc (ASE Calculator object): The original calculator
        outputs (bool): The outputs from the function (Is the calc converged)
        func (str): Path to function that performs the MD like operation
        func_fw_out (str): Path to this function
        func_kwargs (dict): keyword arguments for func
        fw_settings (dict): FireWorks specific settings
    Returns (FWAction): The correct action (restart or updated spec) if convergence is reached
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
        return FWAction(
            update_spec={
                fw_settings["out_spec_atoms"]: atoms,
                fw_settings["out_spec_calc"]: calc,
            }
        )
    del calc["results"]["forces"]
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
    A function that checks if a relaxation is converged (if outputs is True) and either
    stores the relaxed structure in the MongoDB or appends another Firework as its child
    to restart the relaxation
    Args:
        atoms (ASE Atoms object): The original atoms at the start of this job
        calc (ASE Calculator object): The original calculator
        outputs (ASE Atoms Object): The geometry of the final relaxation step
        func (str): Path to function that performs the MD like operation
        func_fw_out (str): Path to this function
        func_kwargs (dict): keyword arguments for func
        fw_settings (dict): FireWorks specific settings
    Returns (FWAction): The correct action (restart or updated spec) if convergence is reached
    """
    func_fw_kwargs["relax_step"] += 1
    aims_out = open(func_kwargs["workdir"] + "/aims.out").readlines()
    converged = "Have a nice day" in aims_out[-2]
    calc = calc2dict(outputs.get_calculator())
    try:
        new_atoms = read_aims(func_kwargs["workdir"] + "/geometry.in.next_step")
        new_atoms.set_calculator(outputs.get_calculator())
        new_atoms_dict = atoms2dict(new_atoms)
        print(new_atoms_dict["cell"])
    except:
        if not converged:
            if "*** WARNING: FHI-aims is terminating due to walltime restrictions\n" in aims_out:
                pass
            else:
                raise IOError(
                    "There was a problem with the FHI Aims calculation stopping program here"
                )
        new_atoms = outputs
        new_atoms_dict = atoms2dict(outputs)
    for key, val in atoms["info"].items():
        if key not in new_atoms_dict["info"]:
            new_atoms_dict["info"][key] = val
    update_spec = dict()
    if converged:
        return FWAction(
            update_spec={
                fw_settings["out_spec_atoms"]: new_atoms_dict,
                fw_settings["out_spec_calc"]: calc,
            }
        )
    else:
        update_spec={
            fw_settings["in_spec_atoms"]: new_atoms_dict,
            fw_settings["in_spec_calc"]: calc,
            fw_settings["kpoint_density_spec"]: k2d(new_atoms, calc["calculator_parameters"]["k_grid"])
        }
    del calc["results"]
    fw_settings["fw_name"] = fw_settings["fw_base_name"] + str(
        func_fw_kwargs["relax_step"]
    )
    fw_settings["spec"].update(update_spec)

    if "to_launchpad" in fw_settings and fw_settings["to_launchpad"]:
        fw_settings["to_launchpad"] = False

    fw = generate_firework(
        func=func,
        func_fw_out=func_fw_out,
        func_kwargs=func_kwargs,
        atoms=new_atoms_dict,
        calc=calc,
        func_fw_out_kwargs=func_fw_kwargs,
        fw_settings=fw_settings,
    )
    return FWAction(detours=[fw], update_spec=update_spec)
