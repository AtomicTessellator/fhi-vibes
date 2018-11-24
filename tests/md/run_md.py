from ase.build import bulk
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.emt import EMT
from ase import units

from hilde.settings import Settings
from hilde.molecular_dynamics import run_md, setup_md


atoms = bulk("Al") * (4, 4, 4)
settings = Settings()

calc = EMT()

MaxwellBoltzmannDistribution(atoms, temp=settings.md.temperature * units.kB)

converged = run_md(atoms, calc, **settings.md)

if not converged:
    from subprocess import run
    run(settings.restart.command.split())
