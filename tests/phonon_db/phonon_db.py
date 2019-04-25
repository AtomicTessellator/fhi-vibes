""" Test for the phonon database """

import numpy as np

from hilde.phonon_db.database_interface import traj_to_database, from_database
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy.postprocess import postprocess as postprocess_ph
from hilde.phono3py.postprocess import postprocess as postprocess_ph3

# get database path
db_path = "test.db"

# Write the initial structure with the phonopy object
hashes_ph = traj_to_database(db_path, "trajectory_phonopy.yaml", True)
# Update the database with third order properties
hashes_ph3 = traj_to_database(db_path, "trajectory_phono3py.yaml", True)

phonon = postprocess_ph(trajectory="trajectory_phonopy.yaml", calculate_full_force_constants=True)
phonon3 = postprocess_ph3(trajectory="trajectory_phono3py.yaml")

ph3_db = from_database(
    db_path,
    get_phonon3=True,
    sc_matrix_2=phonon.get_supercell_matrix(),
    sc_matrix_3=phonon3.get_supercell_matrix(),
    atoms_hash=hashes_ph["atoms_hash"],
    calc_hash=hashes_ph["calc_hash"],
    ph_hash=hashes_ph["ph_hash"],
    ph3_hash=hashes_ph3["ph3_hash"],
)

assert np.max(np.abs(ph3_db.get_fc2()[:] - phonon.get_force_constants()[:])) < 1e-14
assert np.max(np.abs(ph3_db.get_fc3()[:] - phonon3.get_fc3()[:])) < 1e-14

# Get the row from the database
db = connect(db_path)
row = list(
    db.select(
        sc_matrix_2=phonon.get_supercell_matrix(),
        sc_matrix_3=phonon.get_supercell_matrix(),
        atoms_hash=hashes_ph["atoms_hash"],
        calc_hash=hashes_ph["calc_hash"],
        has_fc2=True,
        has_fc3=True,
        columns=["id", "fc_2", "fc_3"],
    )
)[0]

assert np.max(np.abs(row.fc_2[:] - phonon.get_force_constants()[:])) < 1e-14
assert np.max(np.abs(row.fc_3[:] - phonon3.get_fc3()[:])) < 1e-14
