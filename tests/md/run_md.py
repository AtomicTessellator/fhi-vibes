from ase.build import bulk
from ase.calculators.emt import EMT

from hilde import Settings
from hilde.molecular_dynamics import run_md, initialize_md


atoms = bulk("Al") * (4, 4, 4)
settings = Settings()

calc = EMT()

atoms = initialize_md(atoms, **settings.md)

run_md(atoms=atoms, calc=calc, **settings.md)
