from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

from hilde.settings import Settings
from hilde.templates.aims import setup_aims
from hilde.helpers.k_grid import update_k_grid
from hilde.molecular_dynamics import run_md, setup_md


def run(atoms, settings):
    "run and MD"

    calc = setup_aims(settings=settings)

    update_k_grid(atoms, calc, settings.control_kpt.density)

    MaxwellBoltzmannDistribution(atoms, temp=settings.md.temperature * units.kB)

    converged = run_md(atoms, calc, socketio_port=settings.socketio.port, **settings.md)

    return converged


if __name__ == "__main__":
    atoms = read("si.in", format="aims")
    settings = Settings(["../../hilde.cfg", "md.cfg"])

    converged = run(atoms, settings)

    if not converged:
        from subprocess import run

        print("Job not completed, restart.")
        run(settings.restart.command.split())

    print('done.')

