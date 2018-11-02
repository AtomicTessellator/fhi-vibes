""" example on how to use the phonon database """

import numpy as np

from ase.calculators.emt import EMT
from ase.build import bulk

from hilde.helpers.hash import hash_atoms
from hilde.helpers.brillouinzone import get_bands
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy import phono as ph

# Get the settings for the calculation and set up the cell
db_path = "test.db"
print(f"database: {db_path}")

atoms = bulk("Al")
atoms.set_calculator(EMT())

smatrix = 1 * np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
phonon, sc, disp_scs = ph.preprocess(atoms, smatrix)

# connect to the database and check if the calculation was already done
db = connect(db_path)
atoms_hash, calc_hash = hash_atoms(atoms)

force_sets = []
found = False

try:
    rows = list(
        db.select(
            selection=[
                ("supercell_matrix", "=", smatrix),
                ("atoms_hash", "=", atoms_hash),
                ("calc_hash", "=", calc_hash),
                ("is_results", "=", True),
            ]
        )
    )
    if not rows:
        raise KeyError("selection not found")
    else:
        print("Take from database")

except KeyError:
    # if not perform the calculations
    print("not found in database, compute")

    force_sets = [sc.get_forces() for sc in disp_scs]

    phonon.set_forces(force_sets)
    phonon.produce_force_constants()
    phonon.set_band_structure(get_bands(atoms))

    q_mesh = [5, 5, 5]
    phonon.set_mesh(q_mesh)
    phonon.set_total_DOS(freq_pitch=0.1, tetrahedron_method=True)
    phonon.set_thermal_properties(t_step=10, t_max=1000, t_min=0)

    # Write results to the database
    db.write(
        phonon,
        atoms_hash=atoms_hash,
        calc_hash=calc_hash,
        is_results=(phonon.get_force_constants() is not None),
    )

# Example database operations
row = list(
    db.select(
        selection=[
            ("supercell_matrix", "=", smatrix),
            ("atoms_hash", "=", atoms_hash),
            ("calc_hash", "=", calc_hash),
            ("is_results", "=", True),
        ],
        columns=["id", "qmesh", "tp_T", "tp_S", "tp_A", "tp_Cv", "force_constants"],
    )
)[0]

thermalProps = [row.tp_T[30], row.tp_A[30], row.tp_S[30], row.tp_Cv[30]]

force_constants = row.force_constants
phonon = ph.prepare_phonopy(atoms, smatrix, fc2=force_constants)

print(
    f"The thermal properties for this set of calculations "
    + f"(k point mesh {row.qmesh}) are:"
)
print(thermalProps)

# Save some data
force_constants = ph.get_force_constants(phonon)
np.savetxt("force_constants_Al.dat", force_constants)
sc.write("Al.in.supercell")

assert 20 < thermalProps[3] < 30

# print(f"Recalculating the thermal properties with a mesh of [90, 90, 90]")
# phonon = db.get_phonon(
#     selection=[
#         ("supercell_matrix", "=", smatrix),
#         ("atoms_hash", "=", atoms_hash),
#         ("calc_hash", "=", calc_hash),
#         ("is_results", "=", True),
#     ]
# )
# phonon.set_mesh([7, 7, 7])
# phonon.set_thermal_properties(temperatures=row.tp_T)
# print(np.array(phonon.get_thermal_properties()).transpose())
# try:
#     print("The super cell matrix for this calculation was:")
#     print(row.supercell_matrix)
# except:
#     print(
#         f"supercell_matrix was not taken in the initial selection, " + f"resetting rows"
#     )
#     row = list(
#         db.select(
#             selection=[
#                 ("supercell_matrix", "=", smatrix),
#                 ("atoms_hash", "=", atoms_hash),
#                 ("calc_hash", "=", calc_hash),
#                 ("is_results", "=", True),
#             ],
#             columns=["supercell_matrix"],
#         )
#     )[0]
#     print(row.supercell_matrix)
#
