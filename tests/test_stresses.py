import numpy as np
import pytest

from ase.build import bulk
from ase.calculators.lj import LennardJones

from vibes.helpers import get_stresses


@pytest.fixture
def atoms_bulk(stretch=1.5):
    """return a stretched Ar fcc cell"""

    atoms = bulk("Ar", cubic=True)
    atoms.set_cell(atoms.cell * stretch, scale_atoms=True)

    calc = LennardJones(rc=10)
    atoms.calc = calc

    return atoms


def test_get_stresses(atoms_bulk):
    stresses = get_stresses(atoms_bulk)

    assert stresses.shape == (4, 3, 3)

    np.testing.assert_allclose(
        stresses.sum(axis=0), atoms_bulk.get_stress(voigt=False)
    )
