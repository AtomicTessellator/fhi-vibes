""" use the hilde phonopy workflow """

from ase.io import read
from hilde.settings import Settings
from hilde.helpers.k_grid import update_k_grid
from hilde.templates.aims import setup_aims
from hilde.phonopy.workflow import phonopy

atoms = read("si.in")

settings = Settings(["hilde.cfg", "phonopy.cfg"])

calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

phonopy(atoms, calc, socketio_port=settings.socketio.port, **settings.phonopy)
