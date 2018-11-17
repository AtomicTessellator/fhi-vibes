""" use the hilde phonopy workflow """

from ase.io import read

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.helpers.k_grid import update_k_grid
from hilde.phonopy.workflow import phonopy
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import func_to_fire_works
from hilde.templates.aims import setup_aims

atoms = read("si.in")

settings = Settings(["hilde.cfg", "phonopy.cfg"])

calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

if settings.fw_phonon.use_fw:
    fw_settings = settings.fw_phonon
    fw = func_to_fire_works(
        "hilde.phonopy.workflow.phonopy",
        "hilde.tasks.fireworks.fw_action_outs.return_null",
        settings.phonopy,
        atoms,
        calc,
        fw_settings,
        update_calc_settings=None,
    )
    lp = LaunchPadHilde.from_file("/home/purcell/.fireworks/my_launchpad.yaml")
    lp.add_wf(fw)
else:
    phonopy(atoms, calc, socketio_port=settings.socketio.port, **settings.phonopy)
