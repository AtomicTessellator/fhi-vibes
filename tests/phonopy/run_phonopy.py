""" use the hilde phonopy workflow """

from ase.build import bulk
from ase.calculators.emt import EMT
from hilde import Settings
from hilde.phonopy import run_phonopy

atoms = bulk("Al")

calc = EMT()

settings = Settings()

run_phonopy(atoms=atoms, calculator=calc, settings=settings)
