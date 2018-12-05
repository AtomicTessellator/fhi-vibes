from hilde import Settings, setup_aims, restart
from hilde.phonopy.workflow import run


settings = Settings()

atoms, calc = setup_aims(settings=settings)


completed = run(
    atoms, calc, kpt_density=settings.control_kpt.density, **settings.phonopy
)


if not completed:
    restart(settings)
else:
    print("done.")
