"""Functions that generate FWActions after performing Aims Calculations"""
from pathlib import Path

from fireworks import FWAction

from hilde.fireworks.workflow_generator import generate_firework
from hilde.helpers.converters import atoms2dict
from hilde.trajectory import reader as traj_reader


def mod_spec_add(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that appends the current results to a specified spec in the MongoDB
    Args:
        atoms (ASE Atoms): The original atoms at the start of this job
        calc (ASE Calculator): The original calculator
        outputs (dict): The outputs from the function (assumes to be a single bool output)
        func (str): Path to function that performs the MD like operation
        func_fw_out (str): Path to this function
        func_kwargs (dict): keyword arguments for func
        fw_settings (dict): FireWorks specific settings
    Returns (FWAction): Modifies the spec to add the current atoms list to it
    """
    atoms_dict = atoms2dict(outputs)
    return FWAction(mod_spec=[{"_push": {fw_settings["mod_spec_add"]: atoms_dict}}])


def socket_calc_check(func, func_fw_out, *args, fw_settings=None, **kwargs):
    """
    A function that checks if a socket calculation is done, and if not restarts
    Args:
        func (str): Path to function that performs the MD like operation
        func_fw_out (str): Path to this function
        args (list): Arguments passed to the socket calculator function
        fw_settings (dict): FireWorks specific settings
        kwargs (dict): Key word arguments passed to the socket calculator function
    Returns (FWAction): Either a new Firework to restart the calculation or an updated
                        spec with the list of atoms
    """
    update_spec = {
        fw_settings["calc_atoms_spec"]: args[0],
        fw_settings["calc_spec"]: args[1],
        fw_settings["metadata_spec"]: args[2],
    }
    if kwargs["outputs"]:
        if "workdir" in kwargs:
            wd = Path(kwargs["workdir"])
        else:
            wd = Path(".")

        if "trajectory" in kwargs:
            traj = kwargs["trajectory"]
        else:
            traj = "trajectory.yaml"

        ca = traj_reader(str((wd / traj).absolute()), False)
        update_spec[fw_settings["mod_spec_add"]] = [atoms2dict(at) for at in ca]
        return FWAction(update_spec=update_spec)
    fw = generate_firework(
        func="hilde.tasks.calculate.calculate_socket",
        func_fw_out="hilde.tasks.fireworks.fw_out.calculate.socket_calc_check",
        func_kwargs=kwargs,
        atoms_calc_from_spec=False,
        inputs=[
            fw_settings["calc_atoms_spec"],
            fw_settings["calc_spec"],
            fw_settings["metadata_spec"],
        ],
        fw_settings=fw_settings,
    )
    return FWAction(update_spec=update_spec, detours=[fw])
