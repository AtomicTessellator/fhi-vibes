from pathlib import Path
from ase.build import bulk
from vibes.settings import Settings
from vibes.aims.context import AimsContext
from vibes.helpers.hash import hash_atoms, hash_atoms_and_calc

parent = Path(__file__).parent

atoms = bulk("Si") * (2, 2, 2)


def test_hash(atoms=atoms):
    atoms.write("geometry.in", format="aims")

    config_file = parent / "hash.cfg"
    settings = Settings(settings_file=parent / "aims.in", config_file=None)
    settings.machine.basissetloc = parent / settings.machine.basissetloc

    ctx = AimsContext(settings)

    atoms = ctx.ref_atoms
    calc = ctx.get_calculator()

    atoms.calc = calc

    atomshash = hash_atoms(atoms)

    _, calchash = hash_atoms_and_calc(atoms, ignore_file=config_file)

    assert atomshash == "d362270c568a4a9de8a5a867034983c3057c3db0", atomshash
    assert calchash == "edd9178b6838aad8ac6f71f33e49d019a95b0b37", calchash

    Path("geometry.in").unlink()


if __name__ == "__main__":
    test_hash()
