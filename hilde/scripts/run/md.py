from hilde.settings import Settings
from hilde.templates.aims import setup_aims
from hilde.molecular_dynamics import run_md
from hilde.helpers.restarts import restart


settings = Settings()

atoms, calc = setup_aims(settings=settings)


# run the MD
converged = run_md(atoms, calc, **settings.md)


if not converged:
    restart(settings)
else:
    print("done.")
