from hilde import Settings, setup_aims, restart
from hilde.phonopy.workflow import phonopy


settings = Settings()

atoms, calc = setup_aims(settings=settings)


completed = phonopy(
    atoms, calc, kpt_density=settings.control_kpt.density, **settings.phonopy
)


if not completed:
    restart(settings)
else:
    print("done.")
