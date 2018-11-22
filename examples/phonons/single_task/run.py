from subprocess import run

from ase.io import read

from hilde.phonopy.workflow import phonopy
from hilde.settings import Settings, default_config_name
from hilde.templates.aims import setup_aims


def run(atoms, settings):
    "run phonopy"

    calc = setup_aims(settings=settings)

    completed = phonopy(
        atoms,
        calc,
        socketio_port=settings.socketio.port,
        kpt_density=settings.control_kpt.density,
        **settings.phonopy
    )

    return completed


if __name__ == "__main__":

    atoms = read("geometry.in", format="aims")

    settings = Settings(default_config_name)

    converged = run(atoms, settings)

    if not converged:
        if "restart" in settings:
            run(settings.restart.command.split())
        else:
            print("Task not completed, please inspect and rerun.")
    else:
        print('done.')
