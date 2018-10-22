''' An example of how to use FireWorks in conjunction with HilDe'''
from ase.build import bulk
from ase.calculators.emt import EMT

from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import get_smatrix, setup_workdir
from hilde.phonon_db.phonon_db import connect
from hilde.structure.structure import pAtoms, patoms2dict
from hilde.tasks import fireworks as fw
from hilde.helpers.paths import cwd

from fireworks import Firework, LaunchPad, PyTask, Workflow
from fireworks.core.rocket_launcher import rapidfire

db_path = 'test.db'
print(f'database: {db_path}')

atoms = pAtoms(bulk('Ni', 'fcc', a=3.5))
atoms.set_calculator(EMT())

workdir = setup_workdir(atoms)

smatrix = get_smatrix(atoms, n_target=32)

atoms_hash, calc_hash = hash_atoms(atoms)
atoms = patoms2dict(atoms)

db = connect(db_path)
has_fc = False
found = False

# set up the LaunchPad and reset it
launchpad = LaunchPad()
launchpad.reset('', require_password=False)

# create the Firework consisting of a single task
args_init = [atoms, smatrix, workdir]
args_anal = [atoms, smatrix, db.filename]

# Initialize the displacements with phonopy
fw1 = Firework(PyTask({
    "func": fw.initialize_phonopy.name,
    "args": args_init}))

# Calculate the forces for all the displacement cells
fw2 = Firework(PyTask({
    "func": fw.calculate_multiple.name,
    "inputs": ["atom_dicts", "workdirs"]}))

# Calculate the force constants using phonopy
fw3 = Firework(PyTask({
    "func": fw.analyze_phonopy.name,
    "args": args_anal,
    "inputs": ["calc_atoms"]}))

workflow = Workflow([fw1, fw2, fw3], {fw1:[fw2], fw2:[fw3]})

launchpad.add_wf(workflow)

with cwd(workdir / 'fireworks', mkdir=True):
    rapidfire(launchpad)

phonon = db.get_phonon(selection=[("supercell_matrix", "=", smatrix),
                                  ("atoms_hash", "=", atoms_hash),
                                  ("calc_hash", "=", calc_hash),
                                  ("has_fc", "=", True)])

phonon.set_mesh(3 * [3])
_, _, frequencies, _ = phonon.get_mesh()
assert 8 < frequencies.max() < 10

# qpoints, weights, frequencies, _ = phonon.get_mesh()
# for q, w, f in zip(qpoints, weights, frequencies):
#     print(f'q = {q} (weight= {w})')
#     print('# Mode   Frequency')
#     for  ii, fi in enumerate(f):
#         print(f'  {ii+1:3d} {fi:12.7f} THz')
#
# print(f'{frequencies.max()}')
