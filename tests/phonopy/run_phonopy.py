""" use the hilde phonopy workflow """

from ase.build import bulk
from ase.calculators.emt import EMT
from hilde.settings import Settings
from hilde.phonopy.workflow import run

atoms = bulk('Al')

settings = Settings()

calc = EMT()

run(atoms, calc, **settings.phonopy)
