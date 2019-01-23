""" example on how to use the phonon database """

import numpy as np

from hilde.helpers.pickle import pread
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.phonon_db.database_interface import to_database, from_database
from hilde.phonon_db.phonon_db import connect

# Connect to the database
db_path = "test.db"
db = connect(db_path)

# Load the atoms and phonopy objects from the pick file
atoms, phonon = pread("phonopy.pick")

# Get the hashes of the atoms object for easy structure comparison
atoms_hash, calc_hash = hash_atoms_and_calc(atoms)
# Check if the phonopy object is in the database and if not write it to the db

to_database(db_path, obj=phonon, calc=atoms.calc)

ph_db = from_database(db_path, get_phonon=True, sc_matrix_2=phonon.get_supercell_matrix(), atoms_hash=atoms_hash, calc_hash=calc_hash, has_fc2=True)
assert np.max(np.abs(ph_db.get_force_constants() - phonon.get_force_constants()[:])) < 1e-14

# Get the row from the database
row = list(
    db.select(
        sc_matrix_2=phonon.get_supercell_matrix(),
        atoms_hash=atoms_hash,
        calc_hash=calc_hash,
        has_fc2=True,
        columns=["id", "fc_2", "sc_matrix_2", "tp_T", "tp_A"],
    )
)[0]

# Compare the second order force constants and the Helmholtz free energies of the phonopy and row objects
print(
    "Both force constant matrices are the same:",
    np.max(np.abs(row.fc_2[:] - phonon.get_force_constants()[:])) < 1e-14,
)
print(
    "All free energy values are the same:",
    np.all(row.tp_A == phonon.get_thermal_properties()[1]),
)
