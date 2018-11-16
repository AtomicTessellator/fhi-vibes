from pathlib import Path
from ase.build import bulk
from hilde.templates.aims import setup_aims
from hilde.helpers.hash import hash_atoms

atoms = bulk("Si") * (2, 2, 2)

config_file = './tests/hashes/hash.cfg'

calc = setup_aims(config_file=config_file)

atoms.calc = calc

atomshash, calchash = hash_atoms(atoms, ignore_file=config_file)

assert atomshash == "28f353b5fe7a2367f240e1bfa79d8d77b5b3e07e", atomshash
assert calchash == "fe6fd13301cc9523988c9b0df38393cf4fbdd242", calchash
