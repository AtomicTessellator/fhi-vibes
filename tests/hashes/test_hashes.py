from pathlib import Path
from ase.build import bulk
from hilde import Settings
from hilde.aims.context import AimsContext
from hilde.helpers.hash import hash_atoms, hash_atoms_and_calc

parent = Path(__file__).parent

atoms = bulk("Si") * (2, 2, 2)


def test_hash(atoms=atoms):
    atoms.write("geometry.in", format="aims")

    config_file = parent / "hash.cfg"

    ctx = AimsContext(Settings(settings_file=parent / "aims.in", config_file=None))

    atoms = ctx.ref_atoms
    calc = ctx.get_calculator()

    atoms.calc = calc

    atomshash = hash_atoms(atoms)

    _, calchash = hash_atoms_and_calc(atoms, ignore_file=config_file)

    assert atomshash == "b4140109ae1b1fe90f9da3d0ab7e3a76d24b8c8b", atomshash
    assert calchash == "d228b55e927f4f10e15f79feb9b3bf81f997ab7e", calchash

    Path("geometry.in").unlink()


if __name__ == "__main__":
    test_hash()
