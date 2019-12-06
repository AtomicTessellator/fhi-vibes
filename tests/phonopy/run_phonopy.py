""" use the vibes phonopy workflow """

from ase.build import bulk
from ase.calculators.emt import EMT
from vibes import Settings
from vibes.phonopy import run_phonopy
from vibes.phonopy.context import PhonopyContext

atoms = bulk("Al")

calc = EMT()

settings = Settings("phonopy.in")

ctx = PhonopyContext(settings=settings)
ctx.settings.atoms = atoms
ctx.settings.atoms.set_calculator(calc)

run_phonopy(ctx=ctx)
