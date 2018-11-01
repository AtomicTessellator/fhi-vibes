''' An example of how to use FireWorks in conjunction with HilDe'''
from argparse import ArgumentParser
from pathlib import Path

from ase.build import bulk
from ase.calculators.emt import EMT

from fireworks import Firework, LaunchPad, PyTask, Workflow
from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire
from fireworks.utilities.fw_serializers import load_object_from_file

from hilde.fireworks_api_adapter.qlaunch_remote import qlaunch_remote
from hilde.helpers.brillouinzone import get_bands_and_labels
from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import get_smatrix, setup_workdir
from hilde.parsers import read_structure
from hilde.phonon_db.phonon_db import connect
from hilde.settings import Settings
from hilde.structure.structure import pAtoms, patoms2dict
from hilde.tasks import fireworks as fw
from hilde.templates.aims import setup_aims

parser = ArgumentParser()
parser.add_argument("-rh", "--remote_host", nargs="*",
                    help="Remote host to exec qlaunch. Right now, "
                         "only supports running from a config dir.")
parser.add_argument("-ru", "--remote_user",
                    help="Username to login to remote host.")
parser.add_argument("-rp", "--remote_password",
                    help="Password for remote host (if necessary). For "
                         "best operation, it is recommended that you do "
                         "passwordless ssh.")
parser.add_argument("-wd", "--workdir", default=".",
                    help="directory used to calculate the individual atom calculations")
args = parser.parse_args()
# Get the settings for the calculation and set up the cell
st = Settings('../../hilde.cfg')
if st.database.name.startswith('postgres'):
    db_path = st.database.name
else:
    db_path = str(Path(st.database.location) / st.database.name)
print(f'database: {db_path}')

aims_settings = {
    # 'command': st.machine.aims_command,
    'species_dir': str(Path(st.machine.basissetloc) / 'light'),
    'output_level': 'MD_light',
    'relativistic': 'atomic_zora scalar',
    'xc': 'pw-lda',
    'k_grid': 3 * [2],
    "sc_accuracy_rho" : 0.0001,
    "sc_accuracy_forces" : 0.0005
}
aims = setup_aims(aims_settings)
atoms = read_structure('si.in')
atoms.set_calculator(aims)
# workdir = setup_workdir(atoms, "/u/tpurcell/git/hilde/examples/fireworks", False)
workdir = setup_workdir(atoms, args.workdir, False)
smatrix = get_smatrix(atoms, n_target=64)
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
    queueadapter = load_object_from_file(st.FireWorks.queueadapter_file)
    fworker = FWorker.from_file(st.FireWorks.fworker_file) if st.FireWorks.fworker_file else FWorker()

    # create the Firework consisting of a single task
    args_init = [atoms, smatrix, workdir]
    args_fc = [atoms, smatrix]
    args_db = [atoms, db_path]
    # Initialize the displacements with phonopy
    fw1 = Firework(PyTask({"func": fw.initialize_phonopy.name,
                           "args": args_init}))
    fw2 = Firework(PyTask({"func": fw.calculate_multiple.name,
                          "inputs": ["atom_dicts", "workdirs"]}))
    fw3 = Firework(PyTask({"func": fw.calc_phonopy_force_constants.name,
                           "args": args_fc,
                           "inputs": ["calc_atoms"]}))
    fw4 = Firework(PyTask({"func": fw.calc_phonopy_band_structure.name,
                           "args": [],
                           "inputs": ["phonon_dict"]}))
    fw5 = Firework(PyTask({"func": fw.calc_phonopy_dos.name,
                           "args": [[45,45,45]],
                           "inputs": ["phonon_dict"]}))
    fw6 = Firework(PyTask({"func": fw.add_phonon_to_db.name,
                           "args": args_db,
                           "inputs": ["phonon_dict"]}))
    workflow = Workflow([fw1, fw2, fw3, fw4, fw5, fw6], {fw1:[fw2], fw2:[fw3], fw3:[fw4], fw4:[fw5], fw5:[fw6]})
    launchpad.add_wf(workflow)
    rapidfire(launchpad, nlaunches=2)
    qlaunch_remote("rapidfire", maxjobs_queue=250, nlaunches=2, remote_host=args.remote_host,
                   remote_user=args.remote_user, remote_password=args.remote_password,
                   remote_config_dir=["/u/tpurcell/git/hilde/examples/fireworks"], remote_setup=True,
                   reserve=True)
    rapidfire(launchpad, nlaunches=3)

phonon = db.get_phonon(selection=[("supercell_matrix", "=", smatrix),
                                  ("atoms_hash", "=", atoms_hash),
                                  ("calc_hash", "=", calc_hash),
                                  ("has_fc", "=", True)])

bands, labels = get_bands_and_labels(phonon.primitive)
phonon.set_band_structure(bands)
plt = phonon.plot_band_structure(labels=labels)
plt.ylabel('Frequency [THz]')
plt.savefig('phonon_dispersion.pdf')
