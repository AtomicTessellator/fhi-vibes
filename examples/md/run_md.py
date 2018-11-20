from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

from hilde.settings import Settings
from hilde.templates.aims import setup_aims
from hilde.helpers.k_grid import update_k_grid
from hilde.molecular_dynamics import run_md, setup_md


def run(atoms, config_files):
    "run and MD"
    settings = Settings(config_files=config_files)

    calc = setup_aims(settings=settings)

    update_k_grid(atoms, calc, settings.control_kpt.density)

    atoms, md, prepared = setup_md(atoms, **settings.md)

    if not prepared:
        MaxwellBoltzmannDistribution(atoms, temp=settings.md.temperature * units.kB)

    converged = run_md(
        atoms, calc, md, socketio_port=settings.socketio.port, **settings.md
    )

    return converged


if __name__ == "__main__":
    atoms = read("si.in", format="aims")
    config_files = ["hilde.cfg", "md.cfg"]

    run(atoms, config_files)
