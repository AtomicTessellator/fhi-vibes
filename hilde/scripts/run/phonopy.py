import subprocess as sp

from hilde.phonopy.workflow import phonopy
from hilde.settings import Settings
from hilde.templates.aims import setup_aims


settings = Settings()

atoms, calc = setup_aims(settings=settings)


completed = phonopy(
    atoms, calc, kpt_density=settings.control_kpt.density, **settings.phonopy
)


if not completed:
    if "restart" in settings:
        sp.run(settings.restart.command.split())
    else:
        print("Task not completed, please inspect and rerun.")
else:
    print("done.")
