from hilde import Settings, setup_aims, restart
from hilde.phonopy.workflow import run


settings = Settings()

atoms, calc = setup_aims(settings=settings)


completed = run(atoms, calc, **settings.phonopy)


if not completed:
    restart(settings)
else:
    print("done.")
