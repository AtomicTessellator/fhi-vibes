from ase.io import read
from pathlib import Path

from hilde.fireworks.launchpad import LaunchPadHilde
from hilde.helpers.k_grid import update_k_grid
from hilde.relaxation.bfgs import relax
from hilde.settings import Settings
from hilde.fireworks.workflow_generator import generate_firework
from hilde.tasks.fireworks.fw_action_outs import check_relaxation_complete
from hilde.templates.aims import setup_aims

settings = Settings()
fw_settings = dict(settings.fw_settings)
fw_settings["spec"] = {}
atoms, calc = setup_aims(settings=settings)

fw = generate_firework(
    func=relax,
    func_fw_out=check_relaxation_complete,
    func_kwargs=settings.relaxation,
    func_fw_out_kwargs=settings.relaxation,
    atoms=atoms,
    calc=calc,
    fw_settings=fw_settings,
)

lp = LaunchPadHilde.from_file(str(Path.home()/".fireworks/my_launchpad.yaml"))
lp.add_wf(fw)