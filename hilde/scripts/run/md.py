import numpy as np
import subprocess as sp
from pathlib import Path

from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, PhononHarmonics
from ase import units

from hilde.settings import Settings, default_config_name
from hilde.templates.aims import setup_aims
from hilde.molecular_dynamics import run_md, setup_md


settings = Settings(default_config_name)

atoms = read("geometry.in", format="aims")
calc = setup_aims(settings=settings)

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


converged = run_md(atoms, calc, socketio_port=settings.socketio.port, **settings.md)


if not converged:
    if "restart" in settings:
        sp.run(settings.restart.command.split())
    else:
        print("Task not completed, please inspect and rerun.")
else:
    print("done.")
