from ase.symbols import Symbols
from fireworks import FWAction
from pathlib import Path

from hilde.helpers.converters import calc2dict, atoms2dict
from hilde.helpers.fileformats import last_from_yaml
from hilde.phonon_db.database_api import update_phonon_db
from hilde.phonon_db.row import phonon_to_dict
from hilde.tasks.fireworks.general_py_task import generate_firework

mod_name = __name__

def cont_md_out_fw_action(atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings):
    '''
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
    '''
    if "trajectory" in func_kwargs:
        last_step_dict = last_from_yaml(func_kwargs["trajectory"])
    elif "workdir" in func_kwargs:
        last_step_dict = last_from_yaml(Path(func_kwargs["workdir"] + "/trajectory.yaml").absolute())
    else:
        last_step_dict = last_from_yaml(Path("./trajectory.yaml").absolute())

    for key, val in last_step_dict['atoms'].items():
        atoms[key] = val
    calc["results"] = last_step_dict['calculator']
    for key, val in calc.items():
        atoms[key] = val
    next_step = last_step_dict['opt']['nsteps'] + 1

    if outputs:
        if "db_path" in func_fw_kwargs:
            db_path = func_fw_kwargs["db_path"]
            del(func_fw_kwargs["db_path"])
            update_phonon_db(
                db_path,
                atoms,
                atoms,
                **func_fw_kwargs
            )
        return FWAction(update_spec={fw_settings["out_spec_atoms"]: atoms, fw_settings["out_spec_calc"]: calc})
    del(calc['results']['forces'])
    fw_settings["fw_name"] = fw_settings["fw_base_name"]+ str(next_step)
    if "to_launchpad" in fw_settings and fw_settings["to_launcpad"]:
        fw_settings["to_launcpad"] = False
    if "trajectory" in func_kwargs:
        new_traj_list = func_kwargs["trajectory"].split(".")
    elif "workdir" in func_kwargs:
        new_traj_list = str(Path(func_kwargs["workdir"] + "/trajectory.yaml").absolute()).split(".")
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
    fw = generate_firework(func, func_fw_out, func_kwargs, atoms, calc, func_fw_out_kwargs=func_fw_kwargs, fw_settings=fw_settings)
    return FWAction(detours=[fw])

def run_phonopy_fw_out(atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings):
    if outputs:
        return FWAction()
    fw = generate_firework(func, func_fw_out, func_kwargs, atoms, calc, func_fw_out_kwargs=func_fw_kwargs, fw_settings=fw_settings)
    return FWAction(detours=[fw])

def kgrid_opt_fw_out(atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings):
    '''
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
    '''
    last_step_dict = last_from_yaml(func_kwargs["trajectory"])
    for key, val in last_step_dict['atoms'].items():
        atoms[key] = val
    calc["results"] = last_step_dict['calculator']
    for key, val in calc.items():
        atoms[key] = val
    next_step = last_step_dict['opt']['nsteps'] + 1
    if outputs[0]:
        up_spec = {
            fw_settings["out_spec_k_den"]: outputs[1],
            fw_settings["out_spec_atoms"]: atoms,
            fw_settings["out_spec_calc"]: calc2dict(outputs[2]),
        }
        return FWAction(update_spec=up_spec)

    fw_settings["fw_name"] = fw_settings["fw_base_name"] + str(next_step)
    if fw_settings["to_launcpad"]:
        fw_settings["to_launcpad"] = False
    new_traj_list = func_kwargs["trajectory"].split(".")
    try:
        temp_list = new_traj_list[-2].split("_")
        temp_list[-1] = str(int(temp_list[-1]) + 1)
        new_traj_list[-2] = "_".join(temp_list)
        func_kwargs["trajectory"] = ".".join(new_traj_list)
    except:
        new_traj_list[-2] += "_restart_1"
        func_kwargs["trajectory"] = ".".join(new_traj_list)

    fw = generate_firework(func, func_fw_out, func_kwargs, atoms, outputs[2], func_fw_out_kwargs=func_fw_kwargs, fw_settings=fw_settings)
    return FWAction(detours=[fw])

def return_null_atoms(atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings):
    '''
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
    '''
    return FWAction()

def return_null_general(func, func_fw_out, *args, fw_settings=None, **kwargs):
    '''
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
    '''
    return FWAction()

def fw_out_initialize_phonopy(atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings):
    '''
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
    '''
    detours = []
    calc_dict = calc
    if "kpoint_density_spec" in fw_settings:
        del(fw_settings["kpoint_density_spec"])
    for i,sc in enumerate(outputs[2]):
        sc.info["displacement_id"] = i
        sc_dict = atoms2dict(sc)
        for key, val in calc_dict.items():
            sc_dict[key] = val
        calc_kwargs = { 'workdir': func_fw_kwargs['workdir'] + f"/{i:05d}"}
        fw_settings["fw_name"] = f"forces_{Symbols(atoms['numbers']).get_chemical_formula()}_{i}"
        detours.append(
            generate_firework(
                "hilde.tasks.calculate.calculate",
                "hilde.tasks.fireworks.fw_action_outs.mod_spec_add",
                calc_kwargs,
                sc_dict,
                calc_dict,
                atoms_calc_from_spec=False,
                fw_settings=fw_settings,
                update_calc_settings=None,
            )
        )
    return FWAction(update_spec={"phonon": phonon_to_dict(outputs[0])}, detours=detours)

def mod_spec_add(atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings):
    '''
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
    '''
    atoms_dict = atoms2dict(outputs)
    return FWAction(mod_spec=[{"_push": {fw_settings['mod_spec_add']: atoms_dict}}])

def return_up_spec_at_cl_general(func, func_fw_out, func_args, func_kwargs, *args, fw_settings=None, at_spec=None, cl_spec=None, **kwargs):
    '''
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
    '''
    return FWAction(update_spec=up_spec)

