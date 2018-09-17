from ase.calculators.aims import Aims
from pathlib import Path
import numpy as np
from helpers.hash import hash_atoms
from ase.atoms import Atoms
from ase.db import connect
from pprint import pprint
from helpers.paths import cd
from settings import Settings

st = Settings('hilde.conf')
st.print()

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
si = Atoms('Si2',
           positions=[(0.05, 0., 0.),
                      (b/2, b/2, b/2)],
           cell=[[0, b, b],
                 [b, 0, b],
                 [b, b, 0]],
           pbc=True
           )

si.set_calculator(aims)

db = connect(Path(st.database.location) / st.database.name)
atoms_hash, calc_hash = hash_atoms(si, ignore_file='./test/hash_ignore.ini')
print(atoms_hash, calc_hash)

found = False
try:
    si = db.get_atoms(atoms_hash=atoms_hash,
                      calc_hash=calc_hash,
                      is_results=True,
                      attach_calculator=False)
    found = True
except KeyError:
    print('Si not found, will compute')

if input('proceed? ').lower() == 'y':
    pass
else:
    exit()

with cd(aims_tmp_dir):
    tot_en = si.get_total_energy()
print(f'Total energy is {tot_en}')

# # database:
if not found:
    db.write(si,
             atoms_hash=atoms_hash,
             calc_hash=calc_hash,
             is_results=si.calc.results != {})

print('Database:')
for row in db.select():
    print('\n', row['id'])
    pprint(row.__dict__)
