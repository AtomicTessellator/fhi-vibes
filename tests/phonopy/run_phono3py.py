""" use the hilde phonopy workflow """

from ase.build import bulk
from ase.calculators.emt import EMT
from hilde.phono3py import run_phono3py

atoms = bulk("Al")

calc = EMT()

run_phono3py(atoms=atoms, calculator=calc)
