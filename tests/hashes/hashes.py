from pathlib import Path
from ase.build import bulk
from hilde import Settings
from hilde.templates.aims import setup_aims
from hilde.helpers.hash import hash_atoms

atoms = bulk("Si") * (2, 2, 2)

config_file = "./tests/hashes/hash.cfg"

settings = Settings(settings_file=config_file, config_file=None)

calc = setup_aims(settings=settings)

atoms.calc = calc

atomshash, calchash = hash_atoms(atoms, ignore_file=config_file)

assert atomshash == "28f353b5fe7a2367f240e1bfa79d8d77b5b3e07e", atomshash
assert calchash == "da39a3ee5e6b4b0d3255bfef95601890afd80709", calchash
