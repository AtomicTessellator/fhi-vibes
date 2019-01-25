""" example on how to use the phonon database """

import numpy as np

from ase.calculators.emt import EMT
from ase.build import bulk
import importlib as il
from phonopy import Phonopy

from hilde.helpers.hash import hash_atoms_and_calc
from hilde.helpers.brillouinzone import get_bands
from hilde.helpers.supercell import make_cubic_supercell
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy import wrapper as ph
from hilde.tasks.calculate import calculate_multiple
from hilde.phono3py import wrapper as ph3
# from hilde.phono3py.postprocess import get_forces
from hilde.structure.convert import to_Atoms
from hilde.structure.misc import get_sysname

# Get the settings for the calculation and set up the cell
db_path = "test.db"
print(f"database: {db_path}")

atoms = bulk("Al")
atoms.set_calculator(EMT())

_, smatrix2 = make_cubic_supercell(atoms, 256)
_, smatrix3 = make_cubic_supercell(atoms, 32)

q_mesh = [5, 5, 5]

phonopy_settings = {
    'atoms': atoms,
    'supercell_matrix': smatrix2,
    'displacement': 0.01,
    'symprec': 1e-5,
}

phono3py_settings = {
    'atoms': atoms,
    'supercell_matrix': smatrix3,
    'cutoff_pair_distance': 10.0,
    'log_level': 0,
    'displacement': 0.03,
    'q_mesh': q_mesh
}


# connect to the database and check if the calculation was already done
db = connect(db_path)
atoms_hash, calc_hash = hash_atoms_and_calc(atoms)

force_sets = []
found = False

# Check for second order phonons
try:
    rows = list(
        db.select(
            selection=[
                ("sc_matrix_2", "=", smatrix2),
                ("atoms_hash", "=", atoms_hash),
                ("calc_hash", "=", calc_hash),
                ("has_fc2", "=", True),
            ]
        )
    )
    if not rows:
        raise KeyError("selection not found")
    else:
        print("Taking second order force constants from database")
except KeyError:
    print("Atoms with second order force constants not found in the database")
    phonon, sc2, scs2 = ph.preprocess(**phonopy_settings)
    scs2_computed = calculate_multiple(scs2, atoms.calc, f'{get_sysname(atoms)}/fc2')
    fc2_forces = [cell.get_forces() for cell in scs2_computed]
    phonon.produce_force_constants(fc2_forces)
    phonon.set_band_structure(get_bands(atoms))
    phonon.set_mesh(q_mesh)
    phonon.set_thermal_properties()
    phonon.set_total_DOS(freq_pitch=0.1, tetrahedron_method=True)
    # Write Second Order to the Database
    db.write(
        phonon=phonon,
        atoms_hash=atoms_hash,
        calc_hash=calc_hash,
        has_fc2=(phonon.get_force_constants() is not None),
    )

# flokno: Unit Test!
# # Check for third order phonons, while using the previously calculated second order properties
# try:
#     rows = list(
#         db.select(
#             selection=[
#                 ("sc_matrix_3", "=", smatrix3),
#                 ("atoms_hash", "=", atoms_hash),
#                 ("calc_hash", "=", calc_hash),
#                 ("has_fc3", "=", True),
#             ]
#         )
#     )
#     if not rows:
#         raise KeyError("selection not found")
#     else:
#         print("Taking third order force constants from database")
# except KeyError:
#     print("Atoms with third order force_constants not found in database")
#     row = list(
#         db.select(
#             selection=[
#                 ("sc_matrix_2", "=", smatrix2),
#                 ("atoms_hash", "=", atoms_hash),
#                 ("calc_hash", "=", calc_hash),
#                 ("has_fc2", "=", True),
#             ],
#             columns=["id", "fc_2", "sc_matrix_2", "tp_T"],
#         )
#     )[0]
#     phonon3, sc3, scs3 = ph3.preprocess(**phono3py_settings)
#     phonon3.set_fc2(row.fc_2)
# 
#     scs3_computed = calculate_multiple(scs3, atoms.calc, f'{get_sysname(atoms)}/fc3')
#     fc3_forces = get_forces(scs3_computed)
#     phonon3.produce_fc3(fc3_forces)
#     phonon3.run_thermal_conductivity(temperatures=row.tp_T, write_kappa=True)
#     # Update the database with third order properties
#     db.update(
#         row.id,
#         phonon3=phonon3,
#         atoms_hash=atoms_hash,
#         calc_hash=calc_hash,
#         use_second_order=True,
#         has_fc3=(phonon3.get_fc3() is not None),
#     )
# 
# # Example database operations
# row = list(
#     db.select(
#         selection=[
#             ("sc_matrix_2", "=", smatrix2),
#             ("sc_matrix_3", "=", smatrix3),
#             ("atoms_hash", "=", atoms_hash),
#             ("calc_hash", "=", calc_hash),
#             ("has_fc2", "=", True),
#             ("has_fc3", "=", True),
#         ],
#         columns=["id", "qmesh", "tp_T", "tp_S", "tp_A", "tp_Cv", "tp_kappa", "natoms_in_sc_2"],
#     )
# )[0]
# thermalProps = [row.tp_T[30], row.tp_A[30], row.tp_S[30], row.tp_Cv[30], row.tp_kappa[30][0]]
# 
# print(
#     f"The thermal properties for this set of calculations "
#     + f"(k point mesh {row.qmesh}) are:"
# )
# print(thermalProps)
# # Save some data
# phonon3 = db.get_phonon3(selection=[
#             ("sc_matrix_2", "=", smatrix2),
#             ("atoms_hash", "=", atoms_hash),
#             ("calc_hash", "=", calc_hash),
#             ("has_fc2", "=", True),
#         ])
# force_constants = phonon3.get_fc2().swapaxes(1, 2).reshape(2 * (3 * row.natoms_in_sc_2,))
# np.savetxt("force_constants_Al.dat", force_constants)
# to_Atoms(phonon3.get_phonon_supercell()).write("Al.in.supercell_2", format='aims')
# to_Atoms(phonon3.get_supercell()).write("Al.in.supercell_3", format='aims')
# assert 20 < thermalProps[3] < 30
