from ase.build import bulk
from hilde import Settings
from hilde.aims.context import AimsContext
from hilde.helpers.hash import hash_atoms, hash_atoms_and_calc

atoms = bulk("Si") * (2, 2, 2)
atoms.write("geometry.in")

config_file = "hash.cfg"

ctx = AimsContext(Settings(settings_file="aims.in", config_file=None))

atoms = ctx.ref_atoms
calc = ctx.get_calculator()

atoms.calc = calc

atomshash = hash_atoms(atoms)

_, calchash = hash_atoms_and_calc(atoms, ignore_file=config_file)

assert atomshash == "065d7dc8cde297d3691a7b776177780662acdc4c", atomshash
assert calchash == "d228b55e927f4f10e15f79feb9b3bf81f997ab7e", calchash
