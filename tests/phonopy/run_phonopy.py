""" use the hilde phonopy workflow """

from ase.build import bulk
from ase.calculators.emt import EMT
from hilde.phonopy.workflow import run_phonopy

atoms = bulk("Al")

calc = EMT()

run_phonopy(atoms=atoms, calculator=calc)
