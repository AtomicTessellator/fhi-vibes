from ase.build import bulk
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.emt import EMT
from ase import units

from hilde.settings import Settings
from hilde.molecular_dynamics import run_md, setup_md


def run(atoms, config_files):
    "run and MD"
    settings = Settings(config_files=config_files)

    calc = EMT()

    atoms, md, prepared = setup_md(atoms, **settings.md)

    if not prepared:
        MaxwellBoltzmannDistribution(atoms, temp=settings.md.temperature * units.kB)

    converged = run_md(atoms, calc, md, **settings.md)

    return converged


if __name__ == "__main__":
    atoms = bulk("Al") * (4, 4, 4)
    config_files = "md.cfg"

    run(atoms, config_files)
