''' An example of how to use FireWorks in conjunction with HilDe'''
from pathlib import Path

from ase.build import bulk
from ase.calculators.emt import EMT
from ase.db.row import atoms2dict

from fireworks import Firework, LaunchPad, PyTask, Workflow
from fireworks.core.rocket_launcher import rapidfire

from hilde.helpers.hash import hash_atoms
from hilde.helpers.supercell import find_cubic_cell
from hilde.phonon_db.phonon_db import connect
from hilde.settings import Settings
from hilde.structure.structure import pAtoms
from hilde.tasks import fireworks_tasks as fwt

def get_smatrix(at, n_target=64):
    """ Return the supercell matrix for atoms with target size """
    target_size = n_target / len(at)
    return find_cubic_cell(cell=at.cell, target_size=target_size)

def setup_workdir(at):
    """ Set up a working directory """
    vol = at.get_volume()
    wd = Path('./{}/{:.3f}_EMT'.format(at.sysname, vol)).absolute()
    wd.mkdir(parents=True, exist_ok=True)
    return wd

# Get the settings for the calculation and set up the cell
st = Settings('../../hilde.conf')
if st.database.name.startswith('postgres'):
    db_path = st.database.name
else:
    db_path = str(Path(st.database.location) / st.database.name)
print(f'database: {db_path}')

atoms = pAtoms(bulk('Ni', cubic=True))
atoms.set_calculator(EMT())
workdir = setup_workdir(atoms)
smatrix = get_smatrix(atoms, n_target=32)
atoms_hash, calc_hash = hash_atoms(atoms)
atoms = atoms2dict(atoms)

db = connect(db_path)
has_fc = False
found = False
try:
    rows = list(db.select(selection=[("supercell_matrix", "=", smatrix),
                                     ("atoms_hash", "=", atoms_hash),
                                     ("calc_hash", "=", calc_hash)
                                    ]))
    if not rows:
        raise KeyError('selection not found')
    else:
        found = True
    if len(rows) > 1:
        raise ValueError("Multiple rows meet this search criteria, please narrow down your search.")
    if "has_fc" in rows[0] and rows[0].has_fc:
        has_fc = True
except KeyError:
    print("The system was not found in the database")
if not found or not has_fc:
    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset('', require_password=False)

    # create the Firework consisting of a single task
    args_init = [atoms, smatrix, workdir]
    args_anal = [atoms, smatrix, db_path]
    fw1 = Firework(PyTask({"func": fwt.initialize_phonopy_py_task, "args": args_init}))
    fw2 = Firework(PyTask({"func": fwt.analyze_phonopy_py_task,
                           "args":args_anal,
                           "inputs": ["calc_atoms"]}))
    workflow = Workflow([fw1, fw2], {fw1:[fw2]})
    launchpad.add_wf(workflow)
    rapidfire(launchpad)
