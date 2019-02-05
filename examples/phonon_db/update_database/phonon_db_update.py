""" Test for the phonon database """

import numpy as np

from hilde.phonon_db.database_interface import to_database, from_database, traj_to_database
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy.postprocess import postprocess as postprocess_ph
from hilde.phono3py.postprocess import postprocess as postprocess_ph3

# get database path
db_path = "test.db"

# Write the initial structure with the phonopy object
hashes_ph = traj_to_database(db_path, "trajectory_phonopy.yaml", True)
# Update the database with third order properties
hashes_ph3 = traj_to_database(db_path, "trajectory_phono3py.yaml", True)

phonon = postprocess_ph(trajectory="trajectory_phonopy.yaml")
phonon3 = postprocess_ph3(trajectory="trajectory_phono3py.yaml")

ph3_db = from_database(
    db_path,
    get_phonon3=True,
    sc_matrix_2=phonon.get_supercell_matrix(),
    sc_matrix_3=phonon3.get_supercell_matrix(),
    atoms_hash=hashes_ph["atoms_hash"],
    calc_hash=hashes_ph["calc_hash"],
    hashes=[hashes_ph["ph_hash"], hashes_ph3["ph3_hash"]],
)

print(
    "Both second order force constant matrices are the same:",
    np.max(np.abs(ph3_db.get_fc2()[:] - phonon.get_force_constants()[:])) < 1e-14,
)

print(
    "Both third force constant matrices are the same:",
    np.max(np.abs(ph3_db.get_fc3()[:] - phonon3.get_fc3()[:])) < 1e-14,
)