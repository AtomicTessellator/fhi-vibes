""" example on how to use the phonon database """

import numpy as np

from vibes.phonon_db.database_interface import traj_to_database, from_database
from vibes.phonon_db.phonon_db import connect
from vibes.phonopy.postprocess import postprocess

db_path = "test.db"

trajectory_file = "trajectory_phonopy.yaml"

# Add the phonopy calculation contained in the trajectory to the database
phonon_hash = traj_to_database(db_path, trajectory_file)

# Create phonopy object from trajectory directly
phonon = postprocess(trajectory=trajectory_file)

# Retrieve the phonopy object from the database
phonon_db = from_database(db_path, ph_hash=phonon_hash)

# compute difference in force constants to see if they are similar
diff = np.max(np.abs(phonon_db.get_force_constants() - phonon.get_force_constants()[:]))

assert diff < 1e-14

## Manipulate the database directly
# Connect to the database
db = connect(db_path)

# Get the row from the database directly
row = list(db.select(hash=phonon_hash))[0]

# Compare the second order force constants of the phonopy and row objects
diff = np.max(np.abs(row.fc_2[:] - phonon.get_force_constants()[:]))

print("\nBoth force constant matrices are the same:", diff < 1e-14)
