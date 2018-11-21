""" use the hilde phonopy workflow """

from ase.io import read
from hilde.settings import Settings
from hilde.templates.aims import setup_aims
from hilde.phonopy.workflow import phonopy

atoms = read("si.in")

settings = Settings(["hilde.cfg", "phonopy.cfg"])

calc = setup_aims(settings=settings)

phonopy(
    atoms,
    calc,
    socketio_port=settings.socketio.port,
    kpt_density=settings.control_kpt.density,
    backup_settings=settings,
    **settings.phonopy
)
