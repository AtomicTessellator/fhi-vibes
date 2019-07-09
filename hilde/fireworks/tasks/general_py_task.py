"""Standardize python function for FW PyTasks"""
import os

from hilde import DEFAULT_CONFIG_FILE
from hilde.phonon_db.ase_converters import dict2atoms
from hilde.helpers import Timer
from hilde.settings import TaskSettings, Settings


def get_func(func_path):
    """A function that takes in a path to a python function and returns that function

    Parameters
    ----------
    func_path: str
        The path to the python function
    """
    toks = func_path.rsplit(".", 1)
    if len(toks) == 2:
        modname, funcname = toks
        mod = __import__(modname, globals(), locals(), [str(funcname)], 0)
        return getattr(mod, funcname)
    # Handle built in functions.
    return getattr("builtins", toks[0])


def atoms_calculate_task(
    func_path,
    func_fw_out_path,
    func_kwargs,
    func_fw_out_kwargs,
    atoms_dict,
    calc_dict,
    *args,
    fw_settings=None,
):
    """A wrapper function that converts a general function that performs some operation on ASE Atoms/Calculators into a FireWorks style operation

    Parameters
    ----------
    func_path: str
        Path to the function describing the desired set operations to be performed on the Atoms/Calculator objects
    func_fw_out_path: str
        Path to the function that describes how the func inputs/outputs should alter the FireWorks Workflow
    func_kwargs: dict
        A dictionary describing the key word arguments to func
    func_fw_out_kwargs: dict
        Keyword arguments for fw_out function
    atoms_dict: dict
        A dictionary describing the ASE Atoms object
    calc_dict: dict
        A dictionary describing the ASE Calculator
    args: list
        a list of function arguments passed to func
    fw_settings: dict
        A dictionary describing the FireWorks specific settings used in func_fw_out

    Returns
    -------
    FWAction:
        The FWAction func_fw_out outputs

    Raises
    ------
    RuntimeError
        If the Task fails
    """
    start_dir = os.getcwd()
    if fw_settings is None:
        fw_settings = {}

    func = get_func(func_path)
    func_fw_out = get_func(func_fw_out_path)

    default_settings = TaskSettings(name=None, settings=Settings(DEFAULT_CONFIG_FILE))
    calc_dict["command"] = default_settings.machine.aims_command
    if "species_dir" in calc_dict["calculator_parameters"]:
        calc_dict["calculator_parameters"]["species_dir"] = (
            str(default_settings.machine.basissetloc)
            + "/"
            + calc_dict["calculator_parameters"]["species_dir"].split("/")[-1]
        )

    for key, val in calc_dict.items():
        atoms_dict[key] = val
    if "results" in atoms_dict:
        del atoms_dict["results"]
    atoms = dict2atoms(atoms_dict)
    try:
        func_timer = Timer()
        if args:
            outputs = func(atoms, atoms.calc, *args, **func_kwargs)
        else:
            outputs = func(atoms, atoms.calc, **func_kwargs)
        func_fw_out_kwargs["run_time"] = func_timer()
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
    """A wrapper function that converts a general python function into a FireWorks style operation

    Parameters
    ----------
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

    Returns
    -------
    FWAction
        The FWAction func_fw_out outputs
    """
    if fw_settings is None:
        fw_settings = dict()
    func = get_func(func_path)
    func_fw_out = get_func(func_fw_out_path)

    kwargs["outputs"] = func(*args, **kwargs)

    return func_fw_out(
        func_path, func_fw_out_path, *args, fw_settings=fw_settings, **kwargs
    )
