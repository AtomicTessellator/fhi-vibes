#!/usr/bin/env python
# coding: utf-8

# In[1]:


from time import time
from pathlib import Path

from ase.build import bulk
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from ase.calculators.emt import EMT

from hilde.molecular_dynamics import run_md, prepare_from_trajectory

# In[2]:
atoms = bulk("Al") * (4, 4, 4)

# In[3]:
calc = EMT()

# In[4]:
T = 100 * units.kB

md = VelocityVerlet(atoms, timestep=4 * units.fs, logfile="md.log")

trajectory = Path("trajectory.yaml").absolute()

# either take last step or prepare new MD run
if trajectory.exists():
    prepare_from_trajectory(atoms, md)
else:
    MaxwellBoltzmannDistribution(atoms, temp=T)


# In[5]:
max_time = 1

# In[6]:

stime = time()
run_md(atoms, calc, md, walltime=max_time)
print(f"Finished in {time()-stime:.2f}s")

assert time() - stime < max_time