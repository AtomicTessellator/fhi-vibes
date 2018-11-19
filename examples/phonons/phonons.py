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

calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

phonopy(atoms, calc, socketio_port=settings.socketio.port, **settings.phonopy)
