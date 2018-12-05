from hilde import Settings, setup_aims
from hilde.phonopy.workflow import phonopy


settings = Settings()

atoms, calc = setup_aims(settings=settings)


completed = phonopy(atoms, calc, **settings.phonopy)


if not completed:
    print("Task not completed, please inspect and rerun.")
else:
    print("done.")
