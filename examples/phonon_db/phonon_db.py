""" example on how to use the phonon database """

import numpy as np

from hilde.helpers.pickle import pread
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.phonon_db.database_interface import traj_to_database, from_database
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy.postprocess import postprocess
# Connect to the database
db_path = "test.json"
trajectory_file = "trajectory_phonopy.yaml"
# Create phonopy object from trajectory
phonon = postprocess(trajectory=trajectory_file)

# Add the object to the database
ph_hash = traj_to_database(db_path, trajectory_file)
# Retrieve the phonopy object from the database
ph_db = from_database(db_path, ph_hash=ph_hash)

assert np.max(np.abs(ph_db.get_force_constants() - phonon.get_force_constants()[:])) < 1e-14

# Get the row from the database
db = connect(db_path)
row = list(
    db.select(
        sc_matrix_2=phonon.get_supercell_matrix(),
        ph_hash=ph_hash,
        has_fc2=True,
        columns=["id", "fc_2", "sc_matrix_2"],
    )
)[0]

# Compare the second order force constants and the Helmholtz free energies of the phonopy and row objects
print(
    "Both force constant matrices are the same:",
    np.max(np.abs(row.fc_2[:] - phonon.get_force_constants()[:])) < 1e-14,
)