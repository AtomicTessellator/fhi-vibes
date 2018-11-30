from hilde import Settings, setup_aims, restart
from hilde.molecular_dynamics import run_md


settings = Settings()

atoms, calc = setup_aims(settings=settings)


# run the MD
converged = run_md(atoms, calc, **settings.md)


if not converged:
    restart(settings)
else:
    print("done.")
