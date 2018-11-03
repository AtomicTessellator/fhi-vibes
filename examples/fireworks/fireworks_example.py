''' An example of how to use FireWorks in conjunction with HilDe'''
import os
from ase.build import bulk
from ase.calculators.emt import EMT

from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import get_smatrix, setup_workdir
from hilde.phonon_db.phonon_db import connect
from hilde.structure.structure import pAtoms, patoms2dict
from hilde.tasks import fireworks as fw
from hilde.helpers.paths import cwd

from fireworks import Firework, LaunchPad, PyTask, Workflow
from fireworks.core.rocket_launcher import rapidfire, launch_rocket

db_path = (os.getcwd() + '/test.db')
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
try:
    query = {"name": "Ex_WF_Ni", "state": "COMPLETED"}
    wf_ids = launchpad.get_wf_ids(query=query, limit=100)
    for wf_id in wf_ids:
        launchpad.delete_wf(wf_id)
except:
  pass
try:
    query = {"name": "Ex_WF_Ni", "state": "FIZZLED"}
    wf_ids = launchpad.get_wf_ids(query=query, limit=100)
    for wf_id in wf_ids:
        launchpad.delete_wf(wf_id)
except:
    pass
try:
    query = {"name": "Ex_WF_Ni", "state": "READY"}
    wf_ids = launchpad.get_wf_ids(query=query, limit=100)
    for wf_id in wf_ids:
        launchpad.delete_wf(wf_id)
except:
    pass
# create the Firework consisting of a single task
args_init = [smatrix, workdir, atoms]
args_fc = [atoms, smatrix]
args_db = [atoms, db_path]
# Initialize the displacements with phonopy
fw1 = Firework(PyTask({"func": fw.initialize_phonopy.name,
                       "args": args_init}),
               name="initialize_phonopy")
fw2 = Firework(PyTask({"func": fw.calculate_multiple.name,
                       "inputs": ["atom_dicts", "workdirs"]}),
               name="setup_calcs")
fw3 = Firework(PyTask({"func": fw.calc_phonopy_force_constants.name,
                       "args": args_fc,
                       "inputs": ["calc_atoms"]}),
               name="get_fc")
fw4 = Firework(PyTask({"func": fw.add_phonon_to_db.name,
                       "args": args_db,
                       "inputs": ["phonon_dict"]}),
               name="add_to_db")
workflow = Workflow([fw1, fw2, fw3, fw4], {fw1:[fw2], fw2:[fw3], fw3:[fw4]},
                    name="Ex_WF_Ni")
launchpad.add_wf(workflow)
print(workflow)

with cwd(workdir / 'fireworks', mkdir=True):
    launch_rocket(launchpad)
    launch_rocket(launchpad)
    print(workflow)
    launch_rocket(launchpad)
    launch_rocket(launchpad)
    launch_rocket(launchpad)

phonon = db.get_phonon(selection=[("supercell_matrix", "=", smatrix),
                                  ("atoms_hash", "=", atoms_hash),
                                  ("calc_hash", "=", calc_hash),
                                  ("has_fc", "=", True)])

phonon.set_mesh(3 * [3])
_, _, frequencies, _ = phonon.get_mesh()
print(f'Highest frequency: {frequencies.max():.3f} THz (Target: [8,10] THz)')
assert 8 < frequencies.max() < 10
