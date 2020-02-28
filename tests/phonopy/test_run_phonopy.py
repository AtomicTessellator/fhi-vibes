""" use the vibes phonopy workflow """

from pathlib import Path

import pytest

from ase.build import bulk
from ase.calculators.emt import EMT
from vibes import Settings
from vibes.helpers.paths import cwd
from vibes.phono3py.context import Phono3pyContext
from vibes.phonopy.context import PhonopyContext

parent = Path(__file__).parent


atoms = bulk("Al")

calc = EMT()

settings = Settings(parent / "phonopy.in")
ctx = PhonopyContext(settings=settings)

settings3 = Settings(parent / "phono3py.in")
ctx3 = Phono3pyContext(settings=settings3)


@pytest.mark.parametrize("ctx", (ctx, ctx3))
def test_phonopy_ctx(ctx, tmp_path):
    ctx.settings.atoms = atoms
    ctx.settings.atoms.set_calculator(calc)

    with cwd(tmp_path, mkdir=True):
        ctx.run()
