from ase.atoms import Atoms

from fireworks import Firework, PyTask

import os

from hilde import DEFAULT_CONFIG_FILE
from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.helpers.converters import atoms2dict, dict2atoms, calc2dict
from hilde.helpers.k_grid import update_k_grid
from hilde.settings import Settings
from hilde.tasks import fireworks as fw
from hilde.tasks.utility_tasks import update_calc
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

def generate_firework(
    task_spec_list=None,
    atoms=None,
    calc=None,
    atoms_calc_from_spec=False,
    fw_settings=None,
    update_calc_settings={},
    func=None,
    func_fw_out=None,
    func_kwargs=None,
    func_fw_out_kwargs=None,
    args=None,
    inputs=None,
):
    """
    A function that takes in a set of inputs and returns a Firework to perform that operation
    Args:
        task_spec_list (list of TaskSpecs): list of task specifications to perform
        atoms: ASE Atoms object, dictionary or str
            If not atoms__calc_from_spec then this must be an ASE Atoms object or a dictionary describing it
            If atoms__calc_from_spec then this must be a key str to retrieve the Atoms Object from the MongoDB launchpad
        calc: ASE Calculator object, dictionary or str
            If not atoms_calc_from_spec then this must be an ASE Calculator object or a dictionary describing it
            If atoms_calc_from_spec then this must be a key str to retrieve the Calculator from the MongoDB launchpad
        atoms_calc_from_spec: bool
            If True retrieve the atoms/Calculator objects from the MongoDB launchpad
        fw_settings: dict
            Settings used by fireworks to place objects in the right part of the MongoDB
        update_calc_settings: dict
            Used to update the Calculator parameters
        func_fw_out_kwargs: dict
            Keyword functions for the fw_out function
    Returns: Firework
        A Firework that will perform the desired operation on a set of atoms, and process the outputs for Fireworks
    """
    if func:
        if task_spec_list:
            raise ArgumentError("You have defined a task_spec_list and arguments to generate one, please only specify one of these")
        at = atoms is not None
        task_spec_list = [
            TaskSpec(
                func, func_fw_out, at, func_kwargs, func_fw_out_kwargs, args, inputs
            )
        ]
    elif not task_spec_list:
        raise ArgumentError("You have not defined a task_spec_list or arguments to generate one, please specify one of these")

    if "fw_name" not in fw_settings:
        fw_name = None
        fw_settings["fw_base_name"] = ""
    elif "fw_base_name" not in fw_settings:
        fw_settings["fw_base_name"] = fw_settings["fw_name"]

    setup_tasks = []
    if atoms:
        if not atoms_calc_from_spec:
            at = atoms2dict(atoms)
            if "k_grid_density" in update_calc_settings:
                if not isinstance(calc, dict):
                    update_k_grid(atoms, calc, update_calc_settings["k_grid_density"])
                else:
                    recipcell = np.linalg.pinv(at["cell"]).transpose()
                    calc = update_k_grid_calc_dict(calc, recipcell, at["pbc"], new_val)
            cl = calc2dict(calc)
            for key, val in update_calc_settings.items():
                if key != "k_grid_density"
                    cl = update_calc(cl, key, val)
            for key, val in cl.items():
                at[key] = val
        else:
            at = atoms
            cl = calc
            setup_tasks.append(generate_update_calc_task(calc, update_calc_settings))

        if "kpoint_density_spec" in fw_settings:
            setup_tasks.append(generate_mod_calc_task(at, cl, fw_settings["kpoint_density_spec"]))
    else:
        at = None
        cl = None
    job_tasks = []
    for task_spec in task_spec_list:
        job_tasks.append(generate_task(task_spec, fw_settings, at, cl))

    if fw_settings and "to_launcpad" in fw_settings and fw_settings["to_launcpad"]:
        launchpad = LaunchPadHilde.from_file(fw_settings["launchpad_yaml"])
        launchpad.add_wf(Firework(job_tasks, name=fw_settings["fw_name"]))
        return None
    return Firework(setup_tasks + job_tasks, name=fw_settings["fw_name"], spec=fw_settings["spec"])


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
        func,
        func_fw_out,
        at,
        func_kwargs=None,
        func_fw_out_kwargs=None,
        args=None,
        inputs=None,
    ):
        if not isinstance(func, str):
            func = f"{func.__module__}.{func.__name__}"
        if not isinstance(func_fw_out, str):
            func_fw_out = f"{func_fw_out.__module__}.{func_fw_out.__name__}"
        self.func = func
        self.func_fw_out = func_fw_out
        self.at = at

        if func_kwargs:
            self.func_kwargs = func_kwargs
        else:
            self.func_kwargs = {}

        if func_fw_out_kwargs:
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
            self.inputs = {}

    def get_pt_args(self):
        if self.at:
            return [self.func, self.func_fw_out, self.func_kwargs, self.func_fw_out_kwargs]
        return [func, func_fw_out, *args]

    def get_pt_kwargs(self, fw_settings):
        if not fw_settings:
            fw_settings = {}

        if self.at:
            return {"fw_settings": fw_settings}

        to_ret = dict(task_spec.func_kwargs, **task_spec.func_fw_out_kwargs)
        to_ret["fw_settings"] = fw_settings

        return to_ret

    def get_pt_inputs(self):
        if self.at:
            return []
        return inputs