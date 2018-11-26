from hilde.phonopy.workflow import phonopy
from hilde.settings import Settings
from hilde.templates.aims import setup_aims
from hilde.helpers.restarts import restart


settings = Settings()

atoms, calc = setup_aims(settings=settings)


completed = phonopy(
    atoms, calc, kpt_density=settings.control_kpt.density, **settings.phonopy
)


if not completed:
    restart(settings)
else:
    print("done.")
