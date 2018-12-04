""" Example on how to use the ASE (json) database """

from pprint import pprint
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.db import connect
from hilde.helpers.hash import hash_atoms_and_calc


atoms = bulk("Al") * (1, 1, 1)
atoms.calc = EMT()
atoms_hash, calc_hash = hash_atoms_and_calc(atoms)

database_dir = "database.json"
db = connect(database_dir)

found = False
try:
    atoms = db.get_atoms(
        atoms_hash=atoms_hash,
        calc_hash=calc_hash,
        is_results=True,
        attach_calculator=False,
    )
    found = True
except KeyError:
    print("Atoms not found, will compute")

tot_en = atoms.get_total_energy()

print(f"Total energy is {tot_en}")

if not found:
    db.write(
        atoms,
        atoms_hash=atoms_hash,
        calc_hash=calc_hash,
        is_results=atoms.calc.results != {},
    )

print("Database:")
for row in db.select():
    print("\n", row["id"])
    pprint(row.__dict__)
