from ase.io import read

from fireworks import Workflow

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.helpers.k_grid import update_k_grid
from hilde.phonopy.workflow import phonopy
from hilde.relaxation.bfgs import relax
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import atoms_func_to_fireworks, gen_func_to_fireworks
from hilde.templates.aims import setup_aims

atoms = read("geometry.in")

workflow = Settings(["workflow.cfg"])
fw_settings = Settings(['fireworks.cfg'])

fw_list = []

if "relaxation" in workflow:
    relax_settings = Settings(["hilde.cfg"])
    for key, val in  workflow.relaxation_control.items():
        relax_settings.control[key] = val
    calc = setup_aims(settings=relax_settings)
    update_k_grid(atoms, calc, workflow.relaxation_control_kpt.density)
    fw_list.append(
        atoms_func_to_fireworks(
            "hilde.relaxation.bfgs.relax",
            "hilde.tasks.fireworks.fw_action_outs.cont_md_out_fw_action",
            dict(workflow.relaxation),
            atoms,
            calc,
            fw_settings=dict(fw_settings.relaxation_fw),
            update_calc_settings=None,
        )
    )

if "phonopy" in workflow:
    phonopy_settings = Settings(["hilde.cfg"])
    for key, val in workflow.phonopy_control.items():
        phonopy_settings.control[key] = val
    calc = setup_aims(settings=phonopy_settings)
    update_k_grid(atoms, calc, workflow.phonopy_control_kpt.density)
    if fw_settings.phonopy_fw['serial']:
        fw_list.append(
            atoms_func_to_fireworks(
                "hilde.phonopy.workflow.phonopy",
                "hilde.tasks.fireworks.fw_action_outs.return_null_atoms",
                dict(workflow.phonopy),
                fw_settings.phonopy_fw["init_in_spec_atoms"],
                fw_settings.phonopy_fw["init_in_spec_calc"],
                True,
                dict(fw_settings.phonopy_fw),
                update_calc_settings=dict(workflow.phonopy_control),
            )
        )
    else:
        kwargs_init = {"supercell_matrix": workflow.phonopy.supercell_matrix}
        for key in ["displacement", "workdir"]:
            if key in workflow.phonopy:
                kwargs_init[key] = workflow.phonopy[key]

        fw_list.append(
            atoms_func_to_fireworks(
                "hilde.phonopy.workflow.initialize_phonopy_attach_calc",
                "hilde.tasks.fireworks.fw_action_outs.fw_out_initialize_phonopy",
                kwargs_init,
                fw_settings.phonopy_fw["init_in_spec_atoms"],
                fw_settings.phonopy_fw["init_in_spec_calc"],
                atoms_calc_from_spec=True,
                fw_settings=dict(fw_settings.phonopy_fw),
                update_calc_settings=dict(workflow.phonopy_control),
            )
        )
        kwargs = {"fireworks": True}
        anal_keys = ["trajectory", "workdir", "force_constant_file", "displacement"]
        for key in anal_keys:
            if key in workflow.phonopy:
                kwargs[key] = workflow.phonopy[key]
        fw_list.append(
            gen_func_to_fireworks(
                "hilde.phonopy.workflow.analyze_phonopy",
                "hilde.tasks.fireworks.fw_action_outs.return_null_general",
                [],
                ["phonon", fw_settings.phonopy_fw["mod_spec_add"]],
                kwargs,
                dict(fw_settings.phonopy_fw),
            )
        )
fw_dep = {}
for i,fw in enumerate(fw_list[:-1]):
    fw_dep[fw] = fw_list[i+1]
wf = Workflow(fw_list, fw_dep)
lp = LaunchPadHilde(port=27019)
lp.add_wf(wf)