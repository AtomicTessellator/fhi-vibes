from hilde import Settings, setup_aims
from hilde.relaxation.bfgs import relax

settings = Settings()

atoms, calc = setup_aims(settings=settings)

converged = relax(atoms, calc, **settings.relaxation)

if converged:
    print("done.")
else:
    print("Relaxation not converged, please inspect.")
