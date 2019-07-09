"""Creates FireWorks Tasks"""
from fireworks import PyTask


def setup_atoms_task(task_spec, atoms, calc, fw_settings):
    """Setups an ASE Atoms task

    Parameters
    ----------
    task_spec: TaskSpec
        Specification of the Firetask
    atoms: dict
        Dictionary representation of the ase.atoms.Atoms
    calc: dict
        Dictionary representation of the ASE Calculator
    fw_settings: dict
        FireWorks specific parameters

    Returns
    -------
    pt_func: str
        PyTask function name
    pt_args: list
        PyTask args
    pt_inputs: list of str
         PyTask inputs
    pt_kwargs: dict
        PyTask kwargs
    """
    pt_func = "hilde.fireworks.tasks.general_py_task.atoms_calculate_task"
    pt_args = task_spec.pt_args[:4]
    args = task_spec.pt_args[4:]
    pt_inputs = task_spec.pt_inputs
    task_spec.fw_settings = fw_settings
    pt_kwargs = task_spec.pt_kwargs
    if isinstance(atoms, str):
        pt_inputs = [atoms, calc] + pt_inputs
    elif isinstance(calc, str):
        pt_inputs = [calc] + pt_inputs
        pt_args += [atoms]
    else:
        pt_args += [atoms, calc, *args]
    return (pt_func, pt_args, pt_inputs, pt_kwargs)


def setup_general_task(task_spec, fw_settings):
    """Setups a general task

    Parameters
    ----------
    task_spec: TaskSpec
        Specification of the Firetask
    fw_settings: dict
        FireWorks specific parameters

    Returns
    -------
    pt_func: str
        PyTask function name
    pt_args: list
        PyTask args
    pt_inputs: list of str
         PyTask inputs
    pt_kwargs: dict
        PyTask kwargs
    """
    pt_args = task_spec.pt_args
    pt_func = "hilde.fireworks.tasks.general_py_task.general_function_task"
    pt_inputs = task_spec.pt_inputs
    task_spec.fw_settings = fw_settings
    pt_kwargs = task_spec.pt_kwargs
    return (pt_func, pt_args, pt_inputs, pt_kwargs)


def generate_task(task_spec, fw_settings, atoms, calc):
    """Generates a PyTask for a Firework

    Parameters
    ----------
    task_spec: TaskSpec
        Specification of the Firetask
    fw_settings: dict
        FireWorks specific parameters
    atoms: dict
        Dictionary representation of the ase.atoms.Atoms
    calc: dict
        Dictionary representation of the ASE Calculator

    Returns
    -------
    PyTask
        Task for the given TaskSpec
    """
    if task_spec.task_with_atoms_obj:
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
    """Generate a calculator update task

    Parameters
    ----------
    calc_spec: str
        Spec for the calculator in the Fireworks database
    updated_settings: dict
        What parameters to update

    Returns
    -------
    PyTask
        Task to update the calculator in the Fireworks database
    """
    return PyTask(
        {
            "func": "hilde.fireworks.tasks.utility_tasks.update_calc_in_db",
            "args": [calc_spec, updated_settings],
            "inputs": [calc_spec],
        }
    )


def generate_mod_calc_task(at, cl, calc_spec, kpt_spec):
    """Generate a calculator modifier task

    Parameters
    ----------
    at: dict or str
        Either an Atoms dictionary or a spec key to get the Atoms dictionary for the modified system
    cl: dict or str
        Either a Calculator dictionary or a spec key to get the Calculator dictionary for the modified system
    calc_spec: str
        Spec for the calculator in the Fireworks database
    kpt_spec: str
        Spec to update the k-point density of the system

    Returns
    -------
    PyTask
        Task to update the calculator in the Fireworks database
    """
    args = ["k_grid_density", calc_spec]
    kwargs = {"spec_key": kpt_spec}
    if isinstance(cl, str):
        inputs = [cl, kpt_spec]
    else:
        args.append(cl)
        inputs = [kpt_spec]
    if isinstance(at, dict):
        kwargs["atoms"] = at
    else:
        inputs.append(at)
    return PyTask(
        {
            "func": "hilde.fireworks.tasks.utility_tasks.mod_calc",
            "args": args,
            "inputs": inputs,
            "kwargs": kwargs,
        }
    )
