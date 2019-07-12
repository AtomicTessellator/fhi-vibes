""" use the hilde phonopy workflow """

from ase.build import bulk
from ase.calculators.emt import EMT
from hilde import Settings
from hilde.phonopy import run_phonopy
from hilde.phonopy.context import PhonopyContext

atoms = bulk("Al")

calc = EMT()

settings = Settings("phonopy.in")

ctx = PhonopyContext(settings=settings)
ctx.settings.atoms = atoms
ctx.settings.atoms.set_calculator(calc)

run_phonopy(ctx=ctx)
