from fireworks import FWAction

from pathlib import Path

from hilde.fireworks.workflow_generator import generate_firework
from hilde.helpers.converters import calc2dict, atoms2dict, results2dict
from hilde.trajectory import reader as traj_reader
mod_name = __name__


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

def socket_calc_check(
    func, func_fw_out, *args, fw_settings=None, **kwargs
):
    update_spec = {
        fw_settings["calc_atoms_spec"]: args[0],
        fw_settings["calc_spec"]: args[1],
        fw_settings["metadata_spec"]: args[2],
    }
    if kwargs["outputs"]:
        if "workdir" in kwargs:
            wd = Path(kwargs["workdir"])
        else:
            wd = Path('.')

        if "trajectory" in kwargs:
            traj = kwargs["trajectory"]
        else:
            traj = "trajectory.yaml"

        ca = traj_reader(str((wd/traj).absolute()), False)
        update_spec[fw_settings["mod_spec_add"]] = [atoms2dict(at) for at in ca]
        return FWAction(update_spec=update_spec)
    fw = generate_firework(
        func="hilde.tasks.calculate.calculate_socket",
        func_fw_out="hilde.tasks.fireworks.fw_action_outs.calculate.socket_calc_check",
        func_kwargs=kwargs,
        atoms_calc_from_spec=False,
        inputs=[fw_settings["calc_atoms_spec"], fw_settings["calc_spec"], fw_settings["metadata_spec"]],
        fw_settings=fw_settings,
    )
    return FWAction(update_spec=update_spec, detours=[fw])
