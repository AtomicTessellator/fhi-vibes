""" use the hilde phonopy workflow """

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from hilde.statistical_sampling.workflow import run_statistical_sampling

atoms = bulk("Al")
calc = EMT()

run_statistical_sampling(atoms=atoms, calculator=calc)

#r2 = np.genfromtxt("r2.dat")

#assert np.all((r2[:,1] > 0.95) & (r2[:,1] < 1.00))
