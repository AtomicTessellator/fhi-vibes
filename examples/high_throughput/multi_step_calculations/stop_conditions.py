from ase.atoms import Atoms

def test_va(atoms):
    assert isinstance(atoms, Atoms)
    return atoms.get_volume() / len(atoms) < 10
