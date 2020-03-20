""" use the vibes phonopy workflow """

from pathlib import Path

from ase.build import bulk
from ase.calculators.emt import EMT

from vibes import Settings
from vibes.phonopy import run_phonopy
from vibes.phonopy.context import PhonopyContext

parent = Path(__file__).parent


atoms = bulk("Al")

calculator = EMT()

settings = Settings(parent / "phonopy.in")


def test_run_phonopy(atoms=atoms, calc=calculator, settings=settings):
    ctx = PhonopyContext(settings=settings)
    ctx.settings.atoms = atoms
    ctx.settings.atoms.set_calculator(calculator)

    run_phonopy(ctx=ctx)


if __name__ == "__main__":
    test_run_phonopy()
