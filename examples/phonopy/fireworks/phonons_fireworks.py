""" use the hilde phonopy workflow """

from ase.io import read

from fireworks import Workflow

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.helpers.k_grid import update_k_grid
from hilde.phonopy.workflow import phonopy
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import generate_firework
from hilde.templates.aims import setup_aims

settings = Settings()
fw_settings = dict(settings.fw_settings)

atoms, calc = setup_aims(settings=settings)

if "fw_spec" in settings:
    fw_settings["spec"] = dict(settings.fw.fw_spec)
else:
    fw_settings["spec"] = {}
if "fw_spec_qadapter" in settings:
    fw_settings["spec"]["_queueadapter"] = dict(settings.fw_spec_qadapter)

lp = LaunchPadHilde.from_file("/home/purcell/.fireworks/my_launchpad.yaml")
if fw_settings['serial']:
    if "db_storage" in settings  and "db_path" in settings.db_storage:
        settings.phonopy["db_path"] = settings.db_storage.db_path
        settings.phonopy["original_atom_hash"] = atoms_hash
    fw = generate_firework(
        "hilde.phonopy.workflow.phonopy",
        "hilde.tasks.fireworks.fw_action_outs.serial_phonopy_continue",
        dict(settings.phonopy),
        atoms,
        calc,
        atoms_calc_from_spec=False,
        fw_settings=fw_settings,
    )
    lp.add_wf(fw)
else:
    kwargs_init = {"supercell_matrix": settings.phonopy.supercell_matrix}
    kwargs_init_fw_out = {"workdir": settings.phonopy.workdir}
    if "displacement" in settings.phonopy:
        kwargs_init["displacement"] = settings.phonopy.displacement
    init_fw = generate_firework(
        "hilde.phonopy.workflow.initialize_phonopy_attach_calc",
        "hilde.tasks.fireworks.fw_action_outs.add_phonopy_force_calcs",
        kwargs_init,
        atoms,
        calc,
        func_fw_out_kwargs=kwargs_init_fw_out,
        atoms_calc_from_spec=False,
        fw_settings=fw_settings,
    )
    kwargs = {"fireworks": True}
    if "db_storage" in settings  and "db_path" in settings.db_storage:
        kwargs["db_path"] = settings.db_storage.db_path
        kwargs["original_atom_hash"] = atoms_hash

    anal_keys = ["trajectory", "workdir", "force_constant_file", "displacement"]
    for key in anal_keys:
        if key in settings.phonopy:
            kwargs[key] = settings.phonopy[key]
    anal_fw = generate_firework(
        "hilde.phonopy.postprocess.postprocess",
        "hilde.tasks.fireworks.fw_action_outs.fireworks_no_mods_gen_function",
        args=[],
        inputs=["phonon", fw_settings["mod_spec_add"]],
        func_kwargs=kwargs,
        fw_settings=fw_settings,
    )
    wf = Workflow([init_fw, anal_fw], {init_fw: [anal_fw]})
    lp.add_wf(wf)
