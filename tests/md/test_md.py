from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from vibes import Settings
from vibes.molecular_dynamics.context import MDContext
from vibes.konstanten import atomic_units as units

parent = Path(__file__).parent

atoms = bulk("Al") * (4, 4, 4)
settings = Settings(settings_file=parent / "md.in")

calculator = EMT()

np.random.seed(4)
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)


def test_run(tmp_path):
    ctx = MDContext(settings, workdir=tmp_path)

    ctx.atoms = atoms
    ctx.calculator = calculator

    ctx.run()

    # another 5 steps
    ctx.maxsteps += 5
    ctx.run()

    # test log
    assert open(ctx.workdir / "md.log").read() == open(parent / "reference.log").read()


if __name__ == "__main__":
    test_run()
