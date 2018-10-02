from ase.calculators.aims import Aims
from pathlib import Path
import numpy as np
import shutil
from hilde.helpers.hash import hash_atoms
from ase.atoms import Atoms
from ase.io import read, write
from ase.io.jsonio import encode
# from ase.db import connect
from pprint import pprint
from hilde.helpers.paths import cwd
from hilde.settings import Settings
from hilde.parsers import read_aims, read_aims_output
from hilde.phonon_db.phonon_db import connect
from hilde.phonopy import phono as ph
from hilde.structure import pAtoms
from aims_calculation import run_aims


import os

def get_band(q_start, q_end, n_q=51):
    """ Create a list of n_q equidistant qpoints connecting q_start and q_end """
    band = [q_start + (q_end - q_start) / (n_q-1) * i for i in range(n_q)]
    return band


st = Settings('hilde.conf')

aims_out = "aims.out"

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
#

# set up geometry
# Si diamond structure
a = 5.2
b = a/2
siTemp = Atoms('Si2', positions=[(0.1, 0., 0.), (b/2, b/2, b/2)], cell=[[0, b, b], [b, 0, b], [b, b, 0]], pbc=True )
si = pAtoms(ase_atoms=siTemp)
si.set_calculator(aims)

smatrix = 1*np.array([[-1,  1,  1],
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
except KeyError:
    print('Si not found, will compute')

    if input('proceed? ').lower() == 'y':
        pass
    else:
        exit()
    for jj, cell in enumerate(phononCalc[2]):
        folder_with_disp = workdir / f'disp-{jj:03d}'
        folder_with_disp.mkdir(parents=True, exist_ok=True)
        print(str(folder_with_disp) + ' created')
        write(str(folder_with_disp / 'geometry.in'),cell, 'aims', scaled=True)
        try: (folder_with_disp / 'control.in').unlink()
        except: pass
        shutil.copy('control.in', str(folder_with_disp / 'control.in'))
        try:
            forces = read_aims_output(str(folder_with_disp / 'aims.out'))[0].get_forces()
        except FileNotFoundError:
            run_aims(folder_with_disp, "aims.out", st.machine.aims_command)
            forces = read_aims_output(str(folder_with_disp / 'aims.out'))[0].get_forces()
        force_sets.append(forces)
    phonon.set_forces(force_sets)
    phonon.produce_force_constants()
q_mesh = [45, 45, 45]
phonon.set_mesh(q_mesh)
phonon.set_total_DOS(freq_pitch=.1, tetrahedron_method=True )

if not found:
    db.write(phonon, atoms_hash=atoms_hash, calc_hash=calc_hash, is_results=(phonon.get_force_constants() is not None))

# print('Database:')
# for row in db.select():
#     print('\n', row['id'])
#     pprint(row.__dict__)
