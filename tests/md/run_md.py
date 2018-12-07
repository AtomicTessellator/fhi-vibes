from ase.build import bulk
from ase.calculators.emt import EMT

from hilde.settings import Settings
from hilde.molecular_dynamics import run_md, initialize_md
from hilde.helpers.restarts import restart


atoms = bulk("Al") * (4, 4, 4)
settings = Settings()

calc = EMT()

atoms = initialize_md(atoms, **settings.md)

converged = run_md(atoms, calc, **settings.md)

if not converged:
    restart(settings)
