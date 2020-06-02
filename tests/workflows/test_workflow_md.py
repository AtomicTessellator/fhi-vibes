"""test the MD workflow"""

import subprocess as sp
from pathlib import Path

import pytest
from ase.build import bulk

parent = Path(__file__).parent

run_command = "vibes run md"


@pytest.fixture
def atoms():
    """cubic Argon"""
    return bulk("Ar", cubic=True) * (2, 2, 2)


@pytest.mark.parametrize("file", parent.glob("md.*.in"))
def test_npt(atoms, tmp_path, file):
    atoms.write(tmp_path / "geometry.in")
    (tmp_path / "md.in").symlink_to(file)

    sp.run(run_command.split(), cwd=tmp_path)


if __name__ == "__main__":
    test_npt()
