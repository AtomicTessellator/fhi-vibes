from ase.io import read

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.helpers.k_grid import update_k_grid
from hilde.relaxation.bfgs import relax
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import atoms_func_to_fireworks
from hilde.templates.aims import setup_aims

atoms = read("geometry.in")

settings = Settings(["hilde.cfg", "relax.cfg", "fireworks.cfg"])

calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

fw_settings = dict(settings.fw_relax)
fw = atoms_func_to_fireworks(
    "hilde.relaxation.bfgs.relax",
    "hilde.tasks.fireworks.fw_action_outs.cont_md_out_fw_action",
    settings.relaxation,
    atoms,
    calc,
    fw_settings=fw_settings,
    update_calc_settings=None,
)
lp = LaunchPadHilde.from_file("/home/purcell/.fireworks/my_launchpad.yaml")
lp.add_wf(fw)