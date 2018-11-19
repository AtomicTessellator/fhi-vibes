""" use the hilde phonopy workflow """

from ase.io import read

from fireworks import Workflow

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.helpers.k_grid import update_k_grid
from hilde.phonopy.workflow import phonopy
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import atoms_func_to_fireworks, gen_func_to_fireworks
from hilde.templates.aims import setup_aims

atoms = read("si.in")

settings = Settings(["hilde.cfg", "phonopy.cfg"])
fw_settings = dict(Settings(['fireworks.cfg']))['fw_phonon']
calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

lp = LaunchPadHilde.from_file("/home/purcell/.fireworks/my_launchpad.yaml")
if fw_settings['serial']:
    fw = atoms_func_to_fireworks(
        "hilde.phonopy.workflow.phonopy",
        "hilde.tasks.fireworks.fw_action_outs.return_null_atoms",
        settings.phonopy,
        atoms,
        calc,
        fw_settings,
        update_calc_settings=None,
    )
    lp.add_wf(fw)
else:
    init_fw = atoms_func_to_fireworks(
        "hilde.phonopy.workflow.initialize_phonopy",
        "hilde.tasks.fireworks.fw_action_outs.fw_out_initialize_phonopy",
        settings.phonopy,
        atoms,
        calc,
        fw_settings=fw_settings,
        update_calc_settings=None,
    )
    kwargs = settings.phonopy
    kwargs["fireworks"] = True
    anal_fw = gen_func_to_fireworks(
        "hilde.phonopy.workflow.analyze_phonpy",
        "hilde.tasks.fireworks.fw_action_outs.return_null_general",
        [],
        ["phonon", fw_settings["mod_spec_add"]],
        settings.phonopy,
        fw_settings,
    )
    wf = Workflow([init_fw, anal_fw], {init_fw: [anal_fw]})
    lp.add_wf(wf)
