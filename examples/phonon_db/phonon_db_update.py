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


# Load the atoms from the pickle file, and hash the atoms
atoms = pickle.load(open('phonopy.pick', 'rb'))[0]
atoms_hash, calc_hash = hash_atoms_and_calc(atoms)

# connect to the database and check if the calculation was already done
db = connect("phonon.db")

# Pull the Phonopy Object from the database
row = list(
    db.select(
        selection=[
            ("atoms_hash", "=", atoms_hash),
            ("calc_hash", "=", calc_hash),
            ("has_fc2", "=", True),
        ],
        columns=["id", "fc_2", "sc_matrix_2", "tp_T", "tp_A", "tp_S", "tp_Cv"],
    )
)[0]

# Set Phono3py calculation settings

phono3py_settings = {
    'atoms': atoms,
    'fc3_supercell_matrix': [-2, 2, 2, 2, -2, 2, 2, 2, -2],
    'fc2_supercell_matrix': row.sc_matrix_2,
    'log_level': 0
}
# Create a phono3py object and set the second order force constants to those previously calculated
phonon3, _, sc3, _, scs3 = ph3.preprocess(**phono3py_settings)
phonon3.set_fc2(row.fc_2)

# Calculate the forces and produce the third order force constants
scs3_computed = calculate_multiple(scs3, atoms.calc, f'{get_sysname(atoms)}/fc3')
fc3_forces = ph3.get_forces(scs3_computed)
phonon3.produce_fc3(fc3_forces)

# Get the thermal conductivity at the same temperatures as the second order terms were calculated at
phonon3.run_thermal_conductivity(temperatures=row.tp_T, write_kappa=True)
kappa = phonon3.get_thermal_conductivity().get_kappa()[0][30]
header_str = f"{'Temperature':<14}{'A':<10}{'S':<10}{'Cv':<10}{'kappa_xx':<12}{'kappa_yy':<12}{'kappa_zz':<12}{'kappa_yz':<12} {'kappa_xz':<10} kappa_xy"
thermal_prop_str = f"{row.tp_T[30]:<14}{row.tp_A[30]:<10.5}{row.tp_S[30]:<10.5}{row.tp_Cv[30]:<10.5}"
for kk in kappa:
    thermal_prop_str += f"{round(kk,6):<12.5}"
print("The thermal properties at 300 K are:")
print(header_str)
print(thermal_prop_str)
# Update the database with third order properties
db.update(
    row.id,
    phonon3=phonon3,
    atoms_hash=atoms_hash,
    calc_hash=calc_hash,
    use_second_order=True,
    has_fc3=(phonon3.get_fc3() is not None),
)

# Save some data
phonon3_db = db.get_phonon3(
    selection=[
        ("sc_matrix_2", "=", row.sc_matrix_2),
        ("sc_matrix_3", "=", [-2, 2, 2, 2, -2, 2, 2, 2, -2]),
        ("atoms_hash", "=", atoms_hash),
        ("calc_hash", "=", calc_hash),
        ("has_fc2", "=", True),
        ("has_fc3", "=", True),
    ]
)

print("The third order force constant matricies are the same:", np.all(phonon3_db.get_fc3()[:] == phonon3.get_fc3()[:] ) )
kappa_calc = phonon3.get_thermal_conductivity().get_kappa()
kappa_db = phonon3_db.get_thermal_conductivity().get_kappa()
print("The thermal conductivities are the same:", np.all(kappa_calc[:] == kappa_db[:] ) )
