from pathlib import Path
from ase.build import bulk
from hilde.templates.aims import setup_aims
from hilde.helpers.hash import hash_atoms

atoms = bulk("Si") * (2, 2, 2)

calc = setup_aims(config_file="hilde.cfg.template")

atoms.calc = calc

atomshash, calchash = hash_atoms(atoms, ignore_file="./tests/hashes/ignore.cfg")

assert atomshash == "28f353b5fe7a2367f240e1bfa79d8d77b5b3e07e", atomshash
assert calchash == "4016dd5a2d45759f4f8207b62db6e8b917dc4067", calchash
