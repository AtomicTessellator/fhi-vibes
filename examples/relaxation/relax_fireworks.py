from ase.io import read

from hilde.helpers.k_grid import update_k_grid
from hilde.relaxation.bfgs import relax
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import generate_firework
from hilde.tasks.fireworks.fw_action_outs import cont_md_out_fw_action
from hilde.templates.aims import setup_aims

atoms = read("geometry.in")

settings = Settings(["hilde.cfg", "relax.cfg", "fireworks.cfg"])

calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

fw_settings = dict(settings.fw_relax)
generate_firework(
    relax,
    cont_md_out_fw_action,
    settings.relaxation,
    atoms,
    calc,
    fw_settings=fw_settings,
    update_calc_settings=None,
    launchpad="/home/purcell/.fireworks/my_launchpad.yaml"
)
