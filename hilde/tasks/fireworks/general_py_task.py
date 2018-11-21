from ase.atoms import Atoms

from fireworks import Firework, PyTask

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.helpers.converters import atoms2dict, dict2atoms, calc2dict
from hilde.tasks import fireworks as fw

module_name = __name__

def generate_firework(
    func,
    func_fw_out,
    func_kwargs=None,
    func_fw_out_kwargs=None,
    atoms=None,
    calc=None,
    atoms_calc_from_spec=False,
    fw_settings=None,
    update_calc_settings=None,
    args=None,
    inputs=None,
    launchpad=None
):
    '''
    A function that takes in a set of inputs and returns a Firework to perform that operation
    Args:
        func: str or functions
            Path to the function that describes the main operation
            The function must have a setup of (atoms, calc, **kwargs)
        func_fw_out: str or functions
            Path to the function that takes the input and outputs of func and returns the appropriate FWAction for the Workflow
        func_kwargs: dict
            A dictionary that will be passed to func describing its key word arguments
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
    Returns: Firework
        A Firework that will perform the desired operation on a set of atoms, and process the outputs for Fireworks
    '''
    if not isinstance(func, str):
        func = f"{func.__module__}.{func.__name__}"

    if not isinstance(func_fw_out, str):
        func_fw_out = f"{func_fw_out.__module__}.{func_fw_out.__name__}"

    if "fw_name" not in fw_settings:
        fw_name = None
        fw_settings["fw_base_name"] = ""
    elif "fw_base_name" not in fw_settings:
        fw_settings["fw_base_name"] = fw_settings["fw_name"]

    task_list = []
    if atoms:
        if not isinstance(atoms, dict) and not isinstance(atoms, str):
            atoms_dict = atoms2dict(atoms)
        else:
            atoms_dict = atoms

        if not isinstance(calc, dict) and not isinstance(calc, str):
            calc_dict = calc2dict(calc)
        else:
            calc_dict = calc

        if not update_calc_settings:
            update_calc_settings = {}

        pt_args = [
            func,
            func_fw_out,
            func_kwargs,
            func_fw_out_kwargs,
        ]
        if atoms_calc_from_spec:
            pt_inputs = [atoms_dict, calc_dict]
            task_list.append(
                PyTask(
                    {
                        "func": fw.update_calc_in_db.name,
                        "args": [calc_dict, update_calc_settings],
                        "inputs":[calc_dict]
                    }
                )
            )
        else:
            for key, val in update_calc_settings.items():
                if key == "command":
                    calc_dict["command"] = val
                elif key == "basisset_type":
                    sd = calc_dict["calculator_parameters"]["species_dir"].split("/")
                    sd[-1] = val
                    calc_dict["calculator_parameters"]["species_dir"] = "/".join(sd)
                else:
                    if val is None and key in calc_dict["calculator_parameters"]:
                        del(calc_dict["calculator_parameters"][key])
                    elif val is not None:
                        calc_dict["calculator_parameters"][key] = val
            pt_args += [atoms_dict, calc_dict]
            pt_inputs = []
        pt_kwargs = {"fw_settings": fw_settings}
        task_list.append(
            PyTask(
                {
                    "func": fw.atoms_calculate_task.name,
                    "args": pt_args,
                    "inputs": pt_inputs,
                    "kwargs": pt_kwargs
                }
            )
        )
    else:
        kwargs = dict(func_kwargs, **func_fw_out_kwargs)
        if "fw_settings" not in kwargs:
            kwargs["fw_settings"] = fw_settings
        if not args:
            args = []
        task_list.append(
            PyTask(
                {
                    "func": fw.general_function_task.name,
                    "args": [func_path, func_fw_out_path, *args],
                    "inputs": inputs,
                    "kwargs": kwargs,
                }
            )
        )
    if launchpad:
        if isinstance(launchpad, str):
            launchpad = LaunchPadHilde.from_file("/home/purcell/.fireworks/my_launchpad.yaml")
        launchpad.add_wf(Firework(task_list, name=fw_settings["fw_name"]))
        return None
    return Firework(task_list, name=fw_settings["fw_name"])

def get_func(func_path):
    '''A function that takes in a path to a python function and returns that function'''
    toks = func_path.rsplit('.', 1)
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
    fw_settings=None
):
    '''
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
    '''
    if fw_settings is None:
        fw_settings = {}

    func = get_func(func_path)
    func_fw_out = get_func(func_fw_out_path)

    for key, val in calc_dict.items():
        atoms_dict[key] = val
    del(atoms_dict['results'])
    atoms = dict2atoms(atoms_dict)

    outputs = func(atoms, atoms.calc, **func_kwargs)

    return func_fw_out(atoms_dict, calc_dict, outputs, func_path, func_fw_out_path, func_kwargs, func_fw_out_kwargs, fw_settings)

def general_function_task(func_path, func_fw_out_path, *args, fw_settings=None, **kwargs):
    '''
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
    '''
    if fw_settings is None:
        fw_settings = {}
    func = get_func(func_path)
    func_fw_out = get_func(func_fw_out_path)

    outputs = func(*args, **kwargs)

    return func_fw_out(func_path, func_fw_out_path, *args, fw_settings=None, **kwargs)


atoms_calculate_task.name = f"{module_name}.{atoms_calculate_task.__name__}"
general_function_task.name = f"{module_name}.{general_function_task.__name__}"
