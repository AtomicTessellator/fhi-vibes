''' An example of how to use FireWorks in conjunction with HilDe'''
from pathlib import Path

from ase.build import bulk
from ase.calculators.emt import EMT

from fireworks import Firework, LaunchPad, PyTask, Workflow
from fireworks.core.rocket_launcher import rapidfire

from hilde.helpers.brillouinzone import get_bands_and_labels
from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import get_smatrix, setup_workdir
from hilde.phonon_db.phonon_db import connect
from hilde.settings import Settings
from hilde.structure.structure import pAtoms, patoms2dict
from hilde.tasks import fireworks as fw


# Get the settings for the calculation and set up the cell
st = Settings('../../hilde.conf')
if st.database.name.startswith('postgres'):
    db_path = st.database.name
else:
    db_path = str(Path(st.database.location) / st.database.name)
print(f'database: {db_path}')

atoms = pAtoms(bulk('Ni'))
atoms.set_calculator(EMT())
workdir = setup_workdir(atoms)
smatrix = get_smatrix(atoms, n_target=32)
atoms_hash, calc_hash = hash_atoms(atoms)
atoms = patoms2dict(atoms)

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
    # Initialize the displacements with phonopy
    print(fw.initialize_phonopy.name)
    fw1 = Firework(PyTask({"func": fw.initialize_phonopy.name, "args": args_init}))
    # Calculate the forces for all the displacement cells
    fw2 = Firework(PyTask({"func": fw.calculate_multiple.name,
                           "inputs": ["atom_dicts", "workdirs"]}))
    # Calculate the force constants using phonopy
    fw3 = Firework(PyTask({"func": fw.analyze_phonopy.name,
                           "args":args_anal,
                           "inputs": ["calc_atoms"]}))
    workflow = Workflow([fw1, fw2, fw3], {fw1:[fw2], fw2:[fw3]})
    launchpad.add_wf(workflow)
    rapidfire(launchpad)

phonon = db.get_phonon(selection=[("supercell_matrix", "=", smatrix),
                                  ("atoms_hash", "=", atoms_hash),
                                  ("calc_hash", "=", calc_hash),
                                  ("has_fc", "=", True)])

bands, labels = get_bands_and_labels(phonon.primitive)
phonon.set_band_structure(bands)
plt = phonon.plot_band_structure(labels=labels)
plt.ylabel('Frequency [THz]')
plt.savefig('phonon_dispersion.pdf')
