from ase.atoms import Atoms

from fireworks import PyTask

from pathlib import Path

import os

from hilde import DEFAULT_CONFIG_FILE
from hilde.helpers.converters import dict2atoms
from hilde.settings import Settings
from hilde.tasks import fireworks as fw
module_name = __name__

def setup_atoms_task(task_spec, atoms, calc, fw_settings):
    pt_func = fw.atoms_calculate_task.name
    pt_args = task_spec.get_pt_args()
    pt_inputs = task_spec.get_pt_inputs()
    pt_kwargs = task_spec.get_pt_kwargs(fw_settings)
    if isinstance(atoms, str):
        pt_inputs += [atoms, calc]
    else:
        pt_args += [atoms, calc]
    return (pt_func, pt_args, pt_inputs, pt_kwargs)

def setup_general_task(task_spec, fw_settings):
    pt_func = fw.general_function_task.name
    pt_args = task_spec.get_pt_args()
    pt_inputs = task_spec.get_pt_inputs()
    pt_kwargs = task_spec.get_pt_kwargs(fw_settings)
    return (pt_func, pt_args, pt_inputs, pt_kwargs)

def generate_task(task_spec, fw_settings, atoms, calc):
    if task_spec.at:
        pt_params = setup_atoms_task(task_spec, atoms, calc, fw_settings)
    else:
        pt_params = setup_general_task(task_spec, fw_settings)

    return PyTask(
        {
            "func": pt_params[0],
            "args": pt_params[1],
            "inputs": pt_params[2],
            "kwargs": pt_params[3],
        }
    )

def generate_update_calc_task(calc_spec, updated_settings):
    return PyTask(
        {
            "func": fw.update_calc_in_db.name,
            "args": [calc_spec, updated_settings],
            "inputs": [calc_spec],
        }
    )

def generate_mod_calc_task(atoms_dict, calc_spec, kpt_spec):
    return PyTask(
        {
            "func": fw.mod_calc.name,
            "args": ["k_grid_density", calc_spec],
            "inputs": [
                calc_spec,
                kpt_spec,
                atoms_dict,
            ],
            "kwargs": {"spec_key": kpt_spec},
        }
    )

def get_func(func_path):
    """A function that takes in a path to a python function and returns that function"""
    toks = func_path.rsplit(".", 1)
    if len(toks) == 2:
        modname, funcname = toks
        mod = __import__(modname, globals(), locals(), [str(funcname)], 0)
        return getattr(mod, funcname)
    # Handle built in functions.
    return getattr(builtins, toks[0])


def atoms_calculate_task(
    func_path,
    func_fw_out_path,
    func_kwargs,
    func_fw_out_kwargs,
    atoms_dict,
    calc_dict,
    fw_settings=None,
):
    """
    A wrapper function that converts a general function that performs some operation on ASE Atoms/Calculators into a FireWorks style operation
    Args:
        func_path: str
            Path to the function describing the desired set operations to be performed on the Atoms/Calculator objects
        func_fw_out_path: str
            Path to the function that describes how the func inputs/outputs should alter the FireWorks Workflow
        func_kwargs: dict
            A dictionary describing the key word arguments to func
        atoms_dict: dict
            A dictionary describing the ASE Atoms object
        calc_dict: dict
            A dictionary describing the ASE Calculator object
        fw_settings: dict
            A dictionary describing the FireWorks specific settings used in func_fw_out
    Returns: FWAction
        The FWAction func_fw_out outputs
    """
    start_dir = os.getcwd()
    if fw_settings is None:
        fw_settings = {}

    func = get_func(func_path)
    func_fw_out = get_func(func_fw_out_path)

    default_settings = Settings(DEFAULT_CONFIG_FILE)
    if "command" in calc_dict:
        calc_dict["command"] = default_settings.machine.aims_command
    if "species_dir" in calc_dict["calculator_parameters"]:
        calc_dict["calculator_parameters"]["species_dir"] = (
            str(default_settings.machine.basissetloc)
            + "/"
            + calc_dict["calculator_parameters"]["species_dir"].split("/")[-1]
        )

    for key, val in calc_dict.items():
        atoms_dict[key] = val
    del (atoms_dict["results"])
    atoms = dict2atoms(atoms_dict)
    try:
        outputs = func(atoms, atoms.calc, **func_kwargs)
    except:
        os.chdir(start_dir)
        raise RuntimeError(
            f"Function calculation failed, moving to {start_dir} to finish Firework."
        )

    os.chdir(start_dir)
    fw_acts = func_fw_out(
        atoms_dict,
        calc_dict,
        outputs,
        func_path,
        func_fw_out_path,
        func_kwargs,
        func_fw_out_kwargs,
        fw_settings,
    )
    return fw_acts


def general_function_task(
    func_path, func_fw_out_path, *args, fw_settings=None, **kwargs
):
    """
    A wrapper function that converts a general python function into a FireWorks style operation
    Args:
        func_path: str
            Path to the function describing the desired set operations to be performed on the Atoms/Calculator objects
        func_fw_out_path: str
            Path to the function that describes how the func inputs/outputs should alter the FireWorks Workflow
        args: list
            A list of arguments to pass to func and func_fw_out
        fw_settings: dict
            A dictionary describing the FireWorks specific settings used in func_fw_out
        kwargs: dict
            A dict of key word arguments to pass to the func and func_fw_out
    Returns: FWAction
        The FWAction func_fw_out outputs
    """
    if fw_settings is None:
        fw_settings = {}
    func = get_func(func_path)
    func_fw_out = get_func(func_fw_out_path)

    outputs = func(*args, **kwargs)

    return func_fw_out(func_path, func_fw_out_path, *args, fw_settings=None, **kwargs)


atoms_calculate_task.name = f"{module_name}.{atoms_calculate_task.__name__}"
general_function_task.name = f"{module_name}.{general_function_task.__name__}"

class TaskSpec:
    def __init__(
        self,
        func,
        func_fw_out,
        at,
        func_kwargs=None,
        func_fw_out_kwargs=None,
        args=None,
        inputs=None,
        make_abs_path=False
    ):
        if not isinstance(func, str):
            func = f"{func.__module__}.{func.__name__}"
        if not isinstance(func_fw_out, str):
            func_fw_out = f"{func_fw_out.__module__}.{func_fw_out.__name__}"
        self.func = func
        self.func_fw_out = func_fw_out
        self.at = at

        if func_kwargs:
            if "workdir" in func_kwargs and make_abs_path:
                func_kwargs['workdir'] = str(Path(func_kwargs['workdir']).absolute())
            self.func_kwargs = func_kwargs
        else:
            self.func_kwargs = {}

        if func_fw_out_kwargs:
            if "workdir" in func_fw_out_kwargs and make_abs_path:
                func_fw_out_kwargs['workdir'] = str(Path(func_fw_out_kwargs['workdir']).absolute())
            self.func_fw_out_kwargs = func_fw_out_kwargs
        else:
            self.func_fw_out_kwargs = {}

        if args:
            self.args = args
        else:
            self.args = []

        if inputs:
            self.inputs = inputs
        else:
            self.inputs = []

    def get_pt_args(self):
        if self.at:
            return [self.func, self.func_fw_out, self.func_kwargs, self.func_fw_out_kwargs]
        return [self.func, self.func_fw_out, *self.args]

    def get_pt_kwargs(self, fw_settings):
        if not fw_settings:
            fw_settings = {}

        if self.at:
            return {"fw_settings": fw_settings}

        to_ret = dict(self.func_kwargs, **self.func_fw_out_kwargs)
        to_ret["fw_settings"] = fw_settings

        return to_ret

    def get_pt_inputs(self):
        if self.at:
            return []
        return self.inputs