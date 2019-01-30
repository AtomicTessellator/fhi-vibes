from hilde import Settings, setup_aims, restart
from hilde.phono3py.workflow import run


settings = Settings()

atoms, calc = setup_aims(settings=settings)


completed = run(
    atoms, calc, kpt_density=settings.control_kpt.density, **settings.phono3py
)


if not completed:
    restart(settings)
else:
    print("done.")
