from pathlib import Path

from ase.build import bulk

from vibes.calculator.context import CalculatorContext
from vibes.helpers.hash import hash_atoms, hash_atoms_and_calc
from vibes.settings import Settings

parent = Path(__file__).parent

atoms = bulk("Si") * (2, 2, 2)


def test_hash(atoms=atoms):
    atoms.write("geometry.in", format="aims")

    config_file = parent / "hash.cfg"
    settings = Settings(settings_file=parent / "aims.in", config_files=None)
    settings.machine.basissetloc = parent / settings.machine.basissetloc

    ctx = CalculatorContext(settings)

    atoms = ctx.ref_atoms
    calculator = ctx.get_calculator()

    atoms.calc = calculator

    atomshash = hash_atoms(atoms)

    _, calchash = hash_atoms_and_calc(atoms, ignore_file=config_file)

    assert atomshash == "d362270c568a4a9de8a5a867034983c3057c3db0", atomshash
    assert calchash == "6610567789d51f9d79778e4d0304df819b33cb96", calchash

    Path("geometry.in").unlink()


if __name__ == "__main__":
    test_hash()
