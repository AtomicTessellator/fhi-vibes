""" example on how to use the phonon database """

import numpy as np

from hilde.helpers.pickle import pread
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.phonon_db.database_interface import to_database, from_database
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy.postprocess import postprocess
# Connect to the database
db_path = "test.json"
traj = "trajectory_phonopy.yaml"

# Load the atoms and phonopy objects from the pick file
hashes = to_database(db_path, obj=traj)

phonon = postprocess(trajectory=traj)
ph_db = from_database(db_path, identifier=("traj_hash", hashes["traj_hash"]), get_phonon=True)

assert np.max(np.abs(ph_db.get_force_constants() - phonon.get_force_constants()[:])) < 1e-14

# Get the row from the database
db = connect(db_path)
row = list(
    db.select(
        sc_matrix_2=phonon.get_supercell_matrix(),
        atoms_hash=hashes["atoms_hash"],
        calc_hash=hashes["calc_hash"],
        has_fc2=True,
        columns=["id", "fc_2", "sc_matrix_2"],
    )
)[0]

# Compare the second order force constants and the Helmholtz free energies of the phonopy and row objects
print(
    "Both force constant matrices are the same:",
    np.max(np.abs(row.fc_2[:] - phonon.get_force_constants()[:])) < 1e-14,
)