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
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy import wrapper as ph
from hilde.tasks.calculate import calculate_multiple
from hilde.phono3py import wrapper as ph3
from hilde.structure.convert import to_Atoms
from hilde.structure.misc import get_sysname

# Connect to the database
db = connect("test.db")

# Load the atoms and phonopy objects from the pick file
atoms, phonon = pickle.load(open('phonopy.pick', 'rb'))

# Get the hashes of the atoms object for easy structure comparison
atoms_hash, calc_hash = hash_atoms_and_calc(atoms)

# Check if the phonopy object is in the database and if not write it to the db
rows = list(
    db.select(
        selection=[
            ("sc_matrix_2", "=", phonon.get_supercell_matrix()),
            ("atoms_hash", "=", atoms_hash),
            ("calc_hash", "=", calc_hash),
            ("has_fc2", "=", True),
        ],
        columns=["id"],
    )
)
if len(rows) == 0:
    db.write(
        phonon=phonon,
        atoms_hash=atoms_hash,
        calc_hash=calc_hash,
        has_fc2=(phonon.get_force_constants() is not None),
    )

# Get the row from the database
row = list(
    db.select(
        selection=[
            ("sc_matrix_2", "=", phonon.get_supercell_matrix()),
            ("atoms_hash", "=", atoms_hash),
            ("calc_hash", "=", calc_hash),
            ("has_fc2", "=", True),
        ],
        columns=["id", "fc_2", "sc_matrix_2", "tp_T", "tp_A"],
    )
)[0]

# Compare the second order force constants and the Helmholtz free energies of the phonopy and row objects
print("Both force constant matrices are the same:", np.all(row.fc_2[:] == phonon.get_force_constants()[:]))
print("All free energy values are the same:", np.all(row.tp_A == phonon.get_thermal_properties()[1]))
