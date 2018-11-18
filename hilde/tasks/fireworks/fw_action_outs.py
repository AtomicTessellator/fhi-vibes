from fireworks import FWAction

from hilde.helpers.converters import calc2dict, atoms2dict
from hilde.helpers.fileformats import last_from_yaml
from hilde.phonon_db.row import phonon_to_dict
from hilde.tasks.fireworks.general_py_task import atoms_func_to_fireworks
from hilde.tasks.fireworks.general_py_task import gen_func_to_fireworks
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

    fw = atoms_func_to_fireworks(func, func_fw_out, func_kwargs, atoms, calc, fw_settings)
    return FWAction(detours=[fw])

def return_null_atoms(atoms, calc, outputs, func, func_fw_out, func_kwargs, fw_settings):
    return FWAction()

def return_null_general(func, func_fw_out, *args, fw_settings=None, **kwargs):
    return FWAction()

def fw_out_initialize_phonopy(atoms, calc, outputs, func, func_fw_out, func_kwargs, fw_settings):
    detours = []
    calc_dict = calc
    for i,sc in enumerate(outputs[2]):
        sc.info["displacement_id"] = i
        sc_dict = atoms2dict(sc)
        for key, val in calc_dict.items():
            sc_dict[key] = val
        calc_kwargs = { 'workdir': func_kwargs['workdir'] + f"/{i:05d}"}
        detours.append(
            atoms_func_to_fireworks(
                "hilde.tasks.calculate.calculate",
                "hilde.tasks.fireworks.fw_action_outs.mod_spec_add",
                calc_kwargs,
                sc_dict,
                calc_dict,
                fw_settings,
                update_calc_settings=None,
            )
        )
    return FWAction(update_spec={"phonon": phonon_to_dict(outputs[0])}, detours=detours)

def mod_spec_add(atoms, calc, outputs, func, func_fw_out, func_kwargs, fw_settings):
    atoms_dict = atoms2dict(outputs)
    return FWAction(mod_spec=[{"_push": {fw_settings['mod_spec_add']: atoms_dict}}])


