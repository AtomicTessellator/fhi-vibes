""" use the hilde phonopy workflow """

from ase.build import bulk
from ase.calculators.emt import EMT
from hilde.settings import Settings
from hilde.phonopy.workflow import phonopy

atoms = bulk('Al')

settings = Settings(["phonopy.cfg"])

calc = EMT()

phonopy(atoms, calc, **settings.phonopy)
