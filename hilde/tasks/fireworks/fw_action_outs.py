from fireworks import FWAction

from hilde.helpers.fileformats import last_from_yaml
from hilde.tasks.fireworks.general_py_task import func_to_fire_works
mod_name = __name__

def cont_md_out_fw_action(atoms, calc, outputs, func, func_fw_out, func_kwargs, fw_settings):
    last_step_dict = last_from_yaml(func_kwargs["trajectory"])
    for key, val in last_step_dict['atoms'].items():
        atoms[key] = val
    calc["results"] = last_step_dict['calculator']
    for key, val in calc.items():
        atoms[key] = val
    next_step = last_step_dict['opt']['nsteps'] + 1

    if outputs:
        return FWAction(update_spec={fw_settings["out_spec"]: atoms})
    del(calc['results']['forces'])
    fw_settings["fw_name"] = fw_settings["fw_base_name"]+ str(next_step)

    fw = func_to_fire_works(func, func_fw_out, func_kwargs, atoms, calc, fw_settings)
    return FWAction(detours=[fw])

def return_null(atoms, calc, outputs, func, func_fw_out, func_kwargs, fw_settings):
    return FWAction()