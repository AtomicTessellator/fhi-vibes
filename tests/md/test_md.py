import shutil
from pathlib import Path

import numpy as np
from ase import units as u
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from vibes import Settings
from vibes.molecular_dynamics.context import MDContext

parent = Path(__file__).parent

atoms = bulk("Al") * (4, 4, 4)
settings = Settings(settings_file=parent / "md.in")

calc = EMT()

np.random.seed(4)
MaxwellBoltzmannDistribution(atoms, 300 * u.kB)

ctx = MDContext(settings)

ctx.atoms = atoms
ctx.calc = calc


def test_run1():
    ctx.run()


def test_run2():
    # another 5 steps
    ctx.maxsteps += 5
    ctx.run()


def test_log():
    assert open(ctx.workdir / "md.log").read() == open(parent / "reference.log").read()
    shutil.rmtree(ctx.workdir)


if __name__ == "__main__":
    test_run1()
    test_run2()
    test_log()
