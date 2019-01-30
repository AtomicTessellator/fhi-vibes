""" example on how to use the phonon database """

import numpy as np

from ase.calculators.emt import EMT
from ase.build import bulk
import importlib as il
from phonopy import Phonopy

import pickle

from hilde.helpers.hash import hash_atoms_and_calc
from hilde.helpers.brillouinzone import get_bands
from hilde.helpers.supercell import make_cubic_supercell
from hilde.phono3py.postprocess import postprocess as postprocess_ph3
from hilde.phonon_db.database_interface import to_database, from_database
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy import wrapper as ph
from hilde.phonopy.postprocess import postprocess as postprocess_ph
from hilde.tasks.calculate import calculate_multiple
from hilde.phono3py import wrapper as ph3
from hilde.structure.convert import to_Atoms
from hilde.structure.misc import get_sysname

db_path = "test.json"
# Update the database with third order properties
to_database(db_path, "trajectory_phono3py.yaml")

phonon = postprocess_ph(trajectory="trajectory_phonopy.yaml")
phonon3 = postprocess_ph3(trajectory="trajectory_phono3py.yaml")
ph3_db = from_database(
    db_path,
    get_phonon3=True,
    sc_matrix_2=[-1, 1, 1, 1, -1, 1, 1, 1, -1],
    sc_matrix_3=[-1, 1, 1, 1, -1, 1, 1, 1, -1],
)

# print("The atoms_hashes are the same:", (atoms_hash, calc_hash) == hash_atoms_and_calc(at_db))
assert np.max(np.abs(ph3_db.get_fc2()[:] - phonon.get_force_constants()[:]))
assert np.max(np.abs(ph3_db.get_fc3()[:] - phonon3.get_fc3()[:]))
