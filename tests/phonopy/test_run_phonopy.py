"""use the vibes phonopy workflow"""

from pathlib import Path

import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

from vibes import Settings
from vibes.helpers.paths import cwd
from vibes.phonopy.context import PhonopyContext

try:
    from vibes.phono3py.context import Phono3pyContext

    phono3py = True
except ModuleNotFoundError:
    phono3py = False

parent = Path(__file__).parent


atoms = bulk("Al")

calc = EMT()

settings = Settings(parent / "phonopy.in")
ctx = PhonopyContext(settings=settings)

contexts = [ctx]

if phono3py:
    settings3 = Settings(parent / "phono3py.in")
    ctx3 = Phono3pyContext(settings=settings3)
    contexts.append(ctx3)


@pytest.mark.parametrize("ctx", contexts)
def test_phonopy_ctx(ctx, tmp_path):
    ctx.primitive = atoms
    ctx.calculator = calc

    with cwd(tmp_path, mkdir=True):
        ctx.run()
