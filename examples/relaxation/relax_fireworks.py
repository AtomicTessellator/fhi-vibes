from ase.io import read

from hilde.helpers.k_grid import update_k_grid
from hilde.relaxation.bfgs import relax
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import generate_firework
from hilde.tasks.fireworks.fw_action_outs import check_relaxation_complete
from hilde.templates.aims import setup_aims

atoms = read("geometry.in")

settings = Settings(["../../../hilde.cfg", "relax.cfg", "fireworks.cfg"])

calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

fw_settings = dict(settings.fw_relax)
generate_firework(
    relax,
    check_relaxation_complete,
    settings.relaxation,
    atoms,
    calc,
    fw_settings=fw_settings,
    update_calc_settings=None,
)
