from ase import units as u
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.md import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from hilde import Settings
from hilde.molecular_dynamics.workflow import run as run_md


atoms = bulk("Al") * (4, 4, 4)
settings = Settings()

calc = EMT()

MaxwellBoltzmannDistribution(atoms, 300 * u.kB)

md = VelocityVerlet(atoms, timestep=1*u.fs)

print(settings.md)

run_md(atoms=atoms, calc=calc, md=md, **settings.md)
