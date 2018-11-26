import subprocess as sp
from pathlib import Path
import numpy as np

from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, PhononHarmonics
from ase import units

from hilde.settings import Settings
from hilde.templates.aims import setup_aims
from hilde.molecular_dynamics import run_md


settings = Settings()

atoms, calc = setup_aims(settings=settings)


# MD initialization
if not Path("trajectory.yaml").exists():
    if "force_constants" in settings.md:
        print("Initialize positions and velocities using force constants.")
        force_constants = np.loadtxt(settings.md.force_constants)
        PhononHarmonics(
            atoms,
            force_constants,
            quantum=False,
            temp=settings.md.temperature * units.kB,
        )
    else:
        print("Initialize velocities according to Maxwell-Botlzmann distribution.")
        MaxwellBoltzmannDistribution(atoms, temp=settings.md.temperature * units.kB)


# run the MD
converged = run_md(atoms, calc, **settings.md)


if not converged:
    if "restart" in settings:
        sp.run(settings.restart.command.split())
    else:
        print("Task not completed, please inspect and rerun.")
else:
    print("done.")
