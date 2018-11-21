from fireworks import FWAction

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
    last_step_dict = last_from_yaml(func_kwargs["trajectory"])
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

    fw = generate_firework(func, func_fw_out, func_kwargs, func_fw_kwargs, atoms, calc, fw_settings=fw_settings)
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
    for i,sc in enumerate(outputs[2]):
        sc.info["displacement_id"] = i
        sc_dict = atoms2dict(sc)
        for key, val in calc_dict.items():
            sc_dict[key] = val
        calc_kwargs = { 'workdir': func_fw_kwargs['workdir'] + f"/{i:05d}"}
        detours.append(
            generate_firework(
                "hilde.tasks.calculate.calculate",
                "hilde.tasks.fireworks.fw_action_outs.mod_spec_add",
                calc_kwargs,
                calc_kwargs,
                sc_dict,
                calc_dict,
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

