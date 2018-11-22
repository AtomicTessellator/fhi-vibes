from ase.io import read
from hilde.settings import Settings, default_config_name
from hilde.helpers.k_grid import update_k_grid
from hilde.relaxation.bfgs import relax
from hilde.templates.aims import setup_aims

atoms = read("geometry.in")

settings = Settings(default_config_name, write=False)

calc = setup_aims(settings=settings)

update_k_grid(atoms, calc, settings.control_kpt.density)

converged = relax(
    atoms, calc, socketio_port=settings.socketio.port, **settings.relaxation
)

if converged:
    print("done.")
else:
    print("Relaxation not converged, please inspect.")

