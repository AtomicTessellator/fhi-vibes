from ase.symbols import Symbols
from fireworks import FWAction
import numpy as np

from hilde.helpers.converters import calc2dict, atoms2dict
from hilde.fireworks.workflow_generator import generate_firework

mod_name = __name__

def post_bootsrap(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    detours = []
    update_spec = {}
    if "phononpy" in outputs:
        ph_fw_set = fw_settings.copy()
        ph_outputs = outputs["phononpy"]
        ph_settings = func_kwargs["phonon_settings"].copy()
        update_spec["ph_metadata"] = ph_outputs["metadata"]
        ph_fw_set["metadata_spec"] = "ph_metadata"
        calc_dict = calc2dict(ph_outputs["calculator"])
        if ph_settings["serial"]:
            update_spec["ph_calculated_atoms"] = [atoms2dict(at) for at in ph_outputs["atoms_to_calculate"]]
            update_spec["ph_calculator"] = calc_dict
            ph_fw_set["calc_atoms_spec"] = "ph_calculated_atoms"
            ph_fw_set["calc_spec"] = "ph3_calculator"

            detours = add_socket_calc_to_detours(detours, ph_settings, ph_fw_set, "ph")
        else:
            fw_settings["mod_spec_add"] = "ph_forces"
            detours = add_single_calc_to_detours(detours, ph_settings, atoms, ph_outputs["atoms_to_calculate"], calc_dict, fw_settings, "ph")

    if "phonon3py" in outputs:
        ph3_fw_set = fw_settings.copy()
        ph3_outputs = outputs["phonon3py"]
        ph3_settings = func_kwargs["phonon3_settings"].copy()
        update_spec["ph3_metadata"] = ph3_outputs["metadata"]
        ph3_fw_set["metadata_spec"] = "ph3_metadata"
        calc_dict = calc2dict(ph3_outputs["calculator"])
        if ph_settings["serial"]:
            update_spec["ph3_calculated_atoms"] = [atoms2dict(at) for at in ph3_outputs["atoms_to_calculate"]]
            update_spec["ph3_calculator"] = calc_dict
            ph3_fw_set["calc_atoms_spec"] = "ph_calculated_atoms"
            ph3_fw_set["calc_spec"] = "ph3_calculator"
            detours = add_socket_calc_to_detours(detours, ph3_settings, ph3_fw_set, "ph3")
        else:
            fw_settings["mod_spec_add"] = "ph3_forces"
            detours = add_single_calc_to_detours(detours, ph3_settings, atoms, ph3_outputs["atoms_to_calculate"], calc_dict, fw_settings, "ph3")

    return FWAction(update_spec=update_spec, detours=detours)

def add_socket_calc_to_detours(detours, func_kwargs, fw_settings, prefix):
    calc_kwargs = {}
    calc_keys = ["trajectory", "workdir", "backup_folder", "walltime"]
    for key in calc_keys:
        if key in func_fw_kwargs:
            calc_kwargs[key] = func_fw_kwargs[key]
    fw = generate_firework(
        func="hilde.tasks.fireworks.phononpy_phono3py_functions.wrap_calc_socket",
        func_fw_out="hilde.tasks.fireworks.fw_action_outs.calculate.socket_calc_check",
        func_kwargs=calc_kwargs,
        atoms_calc_from_spec=False,
        inputs=[prefix+"_calculated_atoms", prefix+"_calculator", prefix+"_metadata"],
        fw_settings=fw_settings,
    )
    detour.append(fw)

def add_single_calc_to_detours(detours, func_fw_kwargs, atoms, atoms_list, calc_dict, fw_settings, prefix):
    for i, sc in enumerate(atoms_list):
        if not sc:
            continue
        fw_settings=fw_settings.copy()
        sc.info["displacement_id"] = i
        sc_dict = atoms2dict(sc)
        for key, val in calc_dict.items():
            sc_dict[key] = val
        calc_kwargs = {"workdir": func_fw_kwargs["workdir"] + f"/{i:05d}"}
        fw_settings["fw_name"] = prefix + f"forces_{Symbols(atoms['numbers']).get_chemical_formula()}_{i}"
        detours.append(
            generate_firework(
                func="hilde.tasks.calculate.calculate",
                func_fw_out="hilde.tasks.fireworks.fw_action_outs.calculate.mod_spec_add",
                func_kwargs=calc_kwargs,
                atoms=sc_dict,
                calc=calc_dict,
                atoms_calc_from_spec=False,
                fw_settings=fw_settings,
            )
        )
    return detours
