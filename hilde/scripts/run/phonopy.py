import subprocess as sp

from ase.io import read

from hilde.phonopy.workflow import phonopy
from hilde.settings import Settings, default_config_name
from hilde.templates.aims import setup_aims


settings = Settings(default_config_name)

atoms = read("geometry.in", format="aims")
calc = setup_aims(settings=settings)


completed = phonopy(
    atoms,
    calc,
    socketio_port=settings.socketio.port,
    kpt_density=settings.control_kpt.density,
    **settings.phonopy
)


if not completed:
    if "restart" in settings:
        sp.run(settings.restart.command.split())
    else:
        print("Task not completed, please inspect and rerun.")
else:
    print("done.")
