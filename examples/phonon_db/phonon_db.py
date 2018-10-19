from ase.atoms import Atoms
from ase.calculators.aims import Aims
from ase.dft.kpoints import get_cellinfo, special_paths, bandpath

import numpy as np
from pathlib import Path
from pprint import pprint

from hilde.helpers.hash import hash_atoms
from hilde.helpers.brillouinzone import get_bands_and_labels, get_bands
from hilde.parsers.structure import read_structure
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy import phono as ph
from hilde.settings import Settings
from hilde.structure import pAtoms
from hilde.tasks.calculate import calculate_multiple

# Get the settings for the calculation and set up the cell
st = Settings('../../hilde.conf')
if st.database.name.startswith('postgres'):
    db_path = st.database.name
else:
    db_path = str(Path(st.database.location) / st.database.name)
print(f'database: {db_path}')

aims_tmp_dir = Path(st.database.location) / 'aims_tmp'
aims_tmp_dir.mkdir(parents=True, exist_ok=True)

si = read_structure("si.in")
si.set_calculator(Aims(
    aims_command=st.machine.aims_command,
    species_dir=str(Path(st.machine.basissetloc) / 'light'),
    outfilename=str('aims.out'),
    sc_accuracy_rho=1e-4,
    sc_accuracy_forces=5e-4,
    xc='pw-lda',
    k_grid=[2, 2, 2],
    output_level='MD_light'
))

smatrix = 1*np.array([[-1,  1,  1],
                      [ 1, -1,  1],
                      [ 1,  1, -1]])

# connect to the database and check if the calculation was already done
db = connect(db_path)
atoms_hash, calc_hash = hash_atoms(si)
force_sets = []
found = False

try:
    rows = list(db.select(
        selection=[("supercell_matrix", "=", smatrix),
                   ("atoms_hash", "=", atoms_hash),
                   ("calc_hash", "=", calc_hash),
                   ("is_results", "=", True)]))
    if not rows:
        raise KeyError('selection not found')
    else:
        found = True
except KeyError:
    # if not perform the calculations
    print('Si not found, will compute')
    if input('proceed? ').lower() == 'y':
        pass
    else:
        exit()
    phonon, sc, disp_scs = ph.preprocess(si, smatrix)
    vol = si.get_volume()
    workdir = Path(aims_tmp_dir / '{}_{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(
        si.sysname, *smatrix.flatten(), vol))
    # Compute forces and phonon properties
    disp_scs = calculate_multiple(disp_scs, si.calc, workdir, force=True)
    force_sets = [sc.get_forces() for sc in disp_scs]

    phonon.set_forces(force_sets)
    phonon.produce_force_constants()
    phonon.set_band_structure(get_bands(si))

    q_mesh = [45, 45, 45]
    phonon.set_mesh(q_mesh)
    phonon.set_total_DOS(freq_pitch=.1, tetrahedron_method=True )
    phonon.set_thermal_properties(t_step= 25, t_max = 1000, t_min = 0)
    # Write results to the database
    db.write(phonon,
             atoms_hash=atoms_hash,
             calc_hash=calc_hash,
             is_results=(phonon.get_force_constants() is not None))

# Example database operations
row = list(db.select(
    selection=[("supercell_matrix", "=", smatrix),
               ("atoms_hash", "=", atoms_hash),
               ("calc_hash", "=", calc_hash),
               ("is_results", "=", True)],
    columns=["id", "qmesh", "thermal_prop_T", "thermal_prop_S",
             "thermal_prop_A","thermal_prop_Cv" ]))[0]
print(f"The thermal properties for this set of calculations " +
      f"(k point mesh {row.qmesh}) are:")
thermalProps = np.array([row.thermal_prop_T, row.thermal_prop_A,
                         row.thermal_prop_S, row.thermal_prop_Cv]).transpose()
print(thermalProps)
print(f"Recalculating the thermal properties with a mesh of [90, 90, 90]")
phonon = db.get_phonon(
    selection=[("supercell_matrix", "=", smatrix),
               ("atoms_hash", "=", atoms_hash),
               ("calc_hash", "=", calc_hash),
               ("is_results", "=", True)])
phonon.set_mesh([90, 90, 90])
phonon.set_thermal_properties(temperatures=row.thermal_prop_T)
print(np.array(phonon.get_thermal_properties()).transpose())
try:
    print("The super cell matrix for this calculation was:")
    print(row.supercell_matrix)
except:
    print(f"supercell_matrix was not taken in the initial selection, " +
          f"resetting rows")
    row =list(db.select(
    selection=[("supercell_matrix", "=", smatrix),
               ("atoms_hash", "=", atoms_hash),
               ("calc_hash", "=", calc_hash),
               ("is_results", "=", True)],
               columns=["supercell_matrix"]))[0]
    print(row.supercell_matrix)
