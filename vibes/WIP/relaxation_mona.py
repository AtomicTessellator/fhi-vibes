from mona import Rule

from vibes.relaxation.bfgs import relax
from vibes.settings import Settings
from vibes.templates.aims import setup_aims


settings = Settings()

atoms, calc = setup_aims(settings=settings)


@Rule
async def relax_kw(kwargs):
    return relax(**kwargs)


converged = relax(atoms, calc, **settings.relaxation)

if converged:
    print("done.")
else:
    print("Relaxation not converged, please inspect.")