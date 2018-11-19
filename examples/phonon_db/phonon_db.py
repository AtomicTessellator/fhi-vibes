""" example on how to use the phonon database """

import numpy as np

from ase.calculators.emt import EMT
from ase.build import bulk
import importlib as il
from phonopy import Phonopy

from hilde.helpers.hash import hash_atoms
from hilde.helpers.brillouinzone import get_bands
from hilde.helpers.supercell import make_cubic_supercell
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy import phono as ph
from hilde.structure import pAtoms
from hilde.tasks.calculate import calculate_multiple

ph3 = il.import_module('hilde.phono3py.phono3')

# Get the settings for the calculation and set up the cell
db_path = "test.db"
print(f"database: {db_path}")

atoms = pAtoms(bulk("Al"))
atoms.set_calculator(EMT())

_, smatrix2 = make_cubic_supercell(atoms, 32)
_, smatrix3 = make_cubic_supercell(atoms, 32)

q_mesh = [5, 5, 5]
phono3py_settings = {
    'atoms': atoms,
    'fc2_supercell_matrix': smatrix2,
    'fc3_supercell_matrix': smatrix3,
    'cutoff_pair_distance': 3,
    'log_level': 0,
    'q_mesh': q_mesh
}
il.reload(ph3)
phonon3, sc2, sc3, scs2, scs3 = ph3.preprocess(**phono3py_settings)

# connect to the database and check if the calculation was already done
db = connect(db_path)
atoms_hash, calc_hash = hash_atoms(atoms)

force_sets = []
found = False

try:
    rows = list(
        db.select(
            selection=[
                ("sc_matrix_2", "=", smatrix2),
                ("sc_matrix_3", "=", smatrix3),
                ("atoms_hash", "=", atoms_hash),
                ("calc_hash", "=", calc_hash),
                ("has_fc2", "=", True),
                ("has_fc3", "=", True)
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
    scs2_computed = calculate_multiple(scs2, atoms.calc, f'{atoms.sysname}/fc2')
    fc2_forces = ph3.get_forces(scs2_computed)

    scs3_computed = calculate_multiple(scs3, atoms.calc, f'{atoms.sysname}/fc3')
    fc3_forces = ph3.get_forces(scs3_computed)
    print(len(fc2_forces), len(fc3_forces))
    # for force in fc3_forces:
    #     print(force)
    phonon3.produce_fc2(fc2_forces)
    phonon3.produce_fc3(fc3_forces)
    phonon3.run_thermal_conductivity(write_kappa=True)
    phonon = Phonopy(phonon3.get_unitcell(),
                     supercell_matrix=smatrix2,
                     symprec=phonon3._symprec,
                     is_symmetry=phonon3._is_symmetry,
                     factor=phonon3._frequency_factor_to_THz,
                     log_level=phonon3._log_level)
    temps = phonon3.get_thermal_conductivity().get_temperatures()
    phonon.set_force_constants(phonon3.get_fc2())
    phonon.set_band_structure(get_bands(atoms))
    phonon.set_mesh(q_mesh)
    phonon.set_thermal_properties(temperatures=temps)
    phonon.set_total_DOS(freq_pitch=0.1, tetrahedron_method=True)

    # Write results to the database
    db.write(
        phonon3=phonon3,
        phonon=phonon,
        atoms_hash=atoms_hash,
        calc_hash=calc_hash,
        has_fc2=(phonon.get_force_constants() is not None),
        has_fc3=(phonon3.get_fc3() is not None),
    )

# Example database operations
row = list(
    db.select(
        selection=[
            ("sc_matrix_2", "=", smatrix2),
            ("atoms_hash", "=", atoms_hash),
            ("calc_hash", "=", calc_hash),
            ("has_fc2", "=", True),
        ],
        columns=["id", "qmesh", "tp_T", "tp_S", "tp_A", "tp_Cv", "tp_kappa", "fc_2", "fc_3"],
    )
)[0]

thermalProps = [row.tp_T[30], row.tp_A[30], row.tp_S[30], row.tp_Cv[30], row.tp_kappa[30][0]]
force_constants = row.fc_2
# phonon = ph.prepare_phonopy(atoms, smatrix2, fc2=force_constants)

print(
    f"The thermal properties for this set of calculations "
    + f"(k point mesh {row.qmesh}) are:"
)
print(thermalProps)
# Save some data
phonon3 = db.get_phonon3(selection=[
            ("sc_matrix_2", "=", smatrix2),
            ("atoms_hash", "=", atoms_hash),
            ("calc_hash", "=", calc_hash),
            ("has_fc2", "=", True),
        ])
force_constants = phonon3.get_fc2().swapaxes(1, 2).reshape(2 * (3 * row.natoms_in_sc_2,))
np.savetxt("force_constants_Al.dat", force_constants)
sc2.write("Al.in.supercell", format='aims')
assert 20 < thermalProps[3] < 30
