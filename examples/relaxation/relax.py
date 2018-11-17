from ase.io import read

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.helpers.k_grid import update_k_grid
from hilde.relaxation.bfgs import relax
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import func_to_fire_works
from hilde.templates.aims import setup_aims

atoms = read("geometry.in")

settings = Settings(["hilde.cfg", "relax.cfg"])

calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

if not settings.fw_relax.use_fw:
    relax(
        atoms,
        calc,
        socketio_port=settings.socketio.port,
        **settings.relaxation
    )
else:
    fw_settings = settings.fw_relax
    fw = func_to_fire_works(
        "hilde.relaxation.bfgs.relax",
        "hilde.tasks.fireworks.fw_action_outs.cont_md_out_fw_action",
        settings.relaxation,
        atoms,
        calc,
        fw_settings,
        update_calc_settings=None,
    )
    lp = LaunchPadHilde.from_file("/home/purcell/.fireworks/my_launchpad.yaml")
    lp.add_wf(fw)