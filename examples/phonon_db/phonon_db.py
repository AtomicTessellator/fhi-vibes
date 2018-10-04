from ase.calculators.aims import Aims
from pathlib import Path
import numpy as np
from hilde.helpers.hash import hash_atoms
from ase.atoms import Atoms
from ase.dft.kpoints import get_cellinfo, special_paths, bandpath

# from ase.db import connect
from pprint import pprint
from hilde.settings import Settings
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy import phono as ph
from hilde.structure import pAtoms
from hilde.tasks.calculate import compute_forces

st = Settings('hilde.conf')

database_dir = str(Path(st.database.location) / st.database.name)
print(f'database: {database_dir}')
aims_tmp_dir = Path(st.database.location) / 'aims_tmp'
aims_tmp_dir.mkdir(parents=True, exist_ok=True)
aims = Aims(
    aims_command=st.machine.aims_command,
    species_dir=str(Path(st.machine.basissetloc) / 'light'),
    outfilename=str('aims.out'),
    sc_accuracy_rho=1e-4,
    xc='pw-lda',
    k_grid=[2, 2, 2],
    output_level='MD_light'
)

# set up geometry
# Si diamond structure
a = 5.2
b = a/2
siTemp = Atoms('Si2', positions=[(0.1, 0., 0.), (b/2, b/2, b/2)], cell=[[0, b, b], [b, 0, b], [b, b, 0]], pbc=True )
si = pAtoms(ase_atoms=siTemp)
si.set_calculator(aims)

smatrix = 2*np.array([[-1,  1,  1],
                      [ 1, -1,  1],
                      [ 1,  1, -1]])

phononCalc = ph.preprocess(si, smatrix)
phonon = phononCalc[0]
vol = a**3.0

workdir = Path(str(aims_tmp_dir) + '/Si_{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(*smatrix.flatten(), vol))
workdir.mkdir(exist_ok=True)

db = connect(database_dir)
atoms_hash, calc_hash = hash_atoms(si, ignore_file='./test/hash_ignore.ini')

force_sets = []
found = False
try:
    phonon = db.get_phonon(supercell_matrix=smatrix, atoms_hash=atoms_hash, calc_hash=calc_hash, is_results=True, attach_calculator=False)
    found = True
    rows = db.get( selection=None, supercell_matrix=smatrix, atoms_hash=atoms_hash, calc_hash=calc_hash, is_results=True)
    print("The Entropy of Si at 300 K is: ", rows.thermal_entropy(300))
except KeyError:
    print('Si not found, will compute')
    if input('proceed? ').lower() == 'y':
        pass
    else:
        exit()
    phonon.set_forces(compute_forces(phononCalc[2], si.calc, workdir))
    phonon.produce_force_constants()

si_info = get_cellinfo(si.cell)
qpaths = special_paths[si_info.lattice].split(",")
try:
    bands = []
    for path in qpaths:
        for ii, _ in enumerate(path[:-1]):
            bands.append(bandpath([si_info.special_points[path[ii]], si_info.special_points[path[ii+1]]], si.cell)[0])
    phonon.set_band_structure(bands)
except:
    print("Please run the bandstructure calculations.")

q_mesh = [45, 45, 45]
phonon.set_mesh(q_mesh)
phonon.set_total_DOS(freq_pitch=.1, tetrahedron_method=True )
phonon.set_thermal_properties(t_step= 25, t_max = 1000, t_min = 0)

if not found:
    db.write(phonon, atoms_hash=atoms_hash, calc_hash=calc_hash, is_results=(phonon.get_force_constants() is not None))
# print('Database:')
# for row in db.select():
#     print('\n', row['id'])
#     pprint(row.__dict__)
