import pytest
from ase.build import bulk
from vibes.io import get_identifier


@pytest.fixture
def atoms_bulk(type="Al"):
    return bulk(type)


def test_identifier(atoms_bulk):
    id = get_identifier(atoms_bulk)
    assert "space_group" in id
    assert "n_formula" in id
    assert "material" in id
