from ase.build import bulk
from ase.calculators.emt import EMT

from hilde.settings import Settings
from hilde.molecular_dynamics import run_md
from hilde.helpers.restarts import restart


atoms = bulk("Al") * (4, 4, 4)
settings = Settings()

calc = EMT()

converged = run_md(atoms, calc, **settings.md)

if not converged:
    restart(settings)
