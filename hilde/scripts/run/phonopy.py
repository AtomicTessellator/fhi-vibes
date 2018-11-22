import subprocess as sp

from ase.io import read

from hilde.phonopy.workflow import phonopy
from hilde.settings import Settings, default_config_name
from hilde.templates.aims import setup_aims


def run(atoms, calc, settings):
    "run phonopy"


    completed = phonopy(
        atoms,
        calc,
        socketio_port=settings.socketio.port,
        kpt_density=settings.control_kpt.density,
        **settings.phonopy
    )

    return completed


if __name__ == "__main__":

    settings = Settings(default_config_name)

    atoms = read("geometry.in", format="aims")
    calc = setup_aims(settings=settings)

    completed = run(atoms, calc, settings)

    if not completed:
        if "restart" in settings:
            sp.run(settings.restart.command.split())
        else:
            print("Task not completed, please inspect and rerun.")
    else:
        print('done.')
