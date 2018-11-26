from ase.io import read
from pathlib import Path

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.helpers.k_grid import update_k_grid
from hilde.relaxation.bfgs import relax
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import generate_firework
from hilde.tasks.fireworks.fw_action_outs import check_relaxation_complete
from hilde.templates.aims import setup_aims

settings = Settings()
fw_settings = dict(settings.fw_settings)
fw_settings["spec"] = {}
atoms, calc = setup_aims(settings=settings)

fw = generate_firework(
    relax,
    check_relaxation_complete,
    settings.relaxation,
    atoms,
    calc,
    fw_settings=fw_settings,
    update_calc_settings=None,
)

lp = LaunchPadHilde.from_file(str(Path.home()/".fireworks/my_launchpad.yaml"))
lp.add_wf(fw)