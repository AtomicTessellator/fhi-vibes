"""test atoms2json and json2atoms"""
import json
from pathlib import Path

from ase import Atoms
from ase.build import bulk
from vibes.helpers.converters import atoms2json, json2atoms

parent = Path(__file__).parent

atoms = bulk("Si") * (2, 2, 2)
file = parent / "atoms.json"


def test_write(atoms=atoms, file=file):
    """write atoms as json"""
    rep = atoms2json(atoms, reduce=False)
    file.write_text(rep)


def test_read(atoms=atoms, file=file):
    """read atoms as json and compare"""
    rep = file.read_text()
    read_atoms = Atoms(**json.loads(rep))

    assert atoms == read_atoms


def test_write_reduced(atoms=atoms, file=file):
    """write atoms as json with reduced symbols and masses"""
    rep = atoms2json(atoms)
    file.write_text(rep)


def test_read_reduced(atoms=atoms, file=file):
    """read atoms as json and compare"""
    rep = file.read_text()
    read_atoms = json2atoms(rep)

    assert atoms == read_atoms


if __name__ == "__main__":
    test_write()
    test_read()

    test_write_reduced()
    test_read_reduced()
