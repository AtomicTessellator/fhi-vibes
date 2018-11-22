from ase.build import bulk
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.emt import EMT
from ase import units

from hilde.settings import Settings
from hilde.molecular_dynamics import run_md, setup_md


def run(atoms, settings):
    "run and MD"

    calc = EMT()

    MaxwellBoltzmannDistribution(atoms, temp=settings.md.temperature * units.kB)

    converged = run_md(atoms, calc, **settings.md)

    return converged


if __name__ == "__main__":

    atoms = bulk("Al") * (4, 4, 4)

    settings = Settings(config_files="md.cfg")

    converged = run(atoms, settings)

    if not converged:
        from subprocess import run

        run(settings.restart.command.split())
