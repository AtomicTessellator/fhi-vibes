#!/usr/bin/env python
# coding: utf-8

# In[1]:
from time import time
from pathlib import Path
import numpy as np

from ase.build import bulk
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from ase.calculators.emt import EMT

from hilde.templates.aims import setup_aims
from hilde.molecular_dynamics import run_md, prepare_from_trajectory


# In[2]:

# atoms = bulk('Si')
atoms = bulk("Al") * (2, 2, 2)

# In[3]:

port = 12345
aims = setup_aims(
    config_file="../../../hilde.cfg",
    custom_settings={
        "compute_forces": True,
        "k_grid": [2, 2, 2],
        "sc_accuracy_rho": 1e-5,
    },
    port=port,
)

# In[4]:

calc = EMT()
port = None

# In[5]:

T = 100 * units.kB

md = VelocityVerlet(atoms, timestep=4 * units.fs, logfile="md.log")

trajectory = Path("trajectory.yaml").absolute()

# either take last step or prepare new MD run
if trajectory.exists():
    prepare_from_trajectory(atoms, md)
else:
    MaxwellBoltzmannDistribution(atoms, temp=T, rng=np.random.RandomState(1))

# In[7]:

walltime = 1

stime = time()
converged = run_md(
    atoms, calc, md, maxsteps=200, walltime=walltime, workdir="aims", socketio_port=None
)
print(f"Finished in {time()-stime:.2f}s")
print(f".. converged: {converged}")
