from pathlib import Path
from ase.build import bulk
from hilde import Settings
from hilde.templates.aims import setup_aims
from hilde.helpers.hash import hash_atoms, hash_atoms_and_calc

atoms = bulk("Si") * (2, 2, 2)

config_file = "./tests/hashes/hash.cfg"

settings = Settings(settings_file=config_file, config_file=None)

calc = setup_aims(settings=settings)

atoms.calc = calc

atomshash = hash_atoms(atoms)

_, calchash = hash_atoms_and_calc(atoms, ignore_file=config_file)

assert atomshash == "065d7dc8cde297d3691a7b776177780662acdc4c", atomshash
assert calchash == "da39a3ee5e6b4b0d3255bfef95601890afd80709", calchash
