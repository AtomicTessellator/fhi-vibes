from ase.io import read
from hilde.settings import Settings
from hilde.helpers.k_grid import update_k_grid
from hilde.relaxation.bfgs import relax
from hilde.templates.aims import setup_aims

atoms = read("geometry.in")

settings = Settings(["hilde.cfg", "relax.cfg"])

calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

relax(
    atoms,
    calc,
    socketio_port=settings.socketio.port,
    **settings.relaxation
)
