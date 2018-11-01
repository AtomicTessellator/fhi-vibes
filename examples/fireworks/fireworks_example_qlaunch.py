''' An example of how to use FireWorks in conjunction with HilDe'''
from argparse import ArgumentParser
import os
from pathlib import Path

from ase.build import bulk
from ase.calculators.emt import EMT

from fireworks import Firework, LaunchPad, PyTask, Workflow
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
parser.add_argument("-rsd", "--remote_species_dir",
                    help="directory of basis set files for aims on the remote machine")
parser.add_argument("-rcmd", "--remote_command",
                    help="command used to run aims on the remote machine")
parser.add_argument("-wd", "--workdir", default=".",
                    help="directory used to calculate the individual atom calculations")
parser.add_argument("--no_kerberos", action="store_true", help="If set do not use gss_api authentication")
args = parser.parse_args()

# Get the settings for the calculation and set up the cell
db_path = os.getcwd() + '/test.db'
print(f'database: {db_path}')

aims_settings = {
    'species_type' : "light",
    'output_level': 'MD_light',
    'relativistic': 'atomic_zora scalar',
    'xc': 'pw-lda',
    'k_grid': 3 * [2],
    "sc_accuracy_rho" : 0.0001,
    "sc_accuracy_forces" : 0.0005
}
aims = setup_aims(aims_settings)
atoms = read_structure('si.in')
# Distortion done to illustrate submitting multiple jobs to the queue at once
atoms.positions[0][0] += 1e-2
atoms.set_calculator(aims)
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
    # create the Firework consisting of a single task
    args_init = [atoms, smatrix, workdir]
    calc_mods = {"command": args.remote_command, "species_dir": args.remote_species_dir}
    spec_qad = {"_queueadapter": {"nodes": 1, "queue": "express", "walltime": "00:05:00"}}
    kwargs_cm = { "calc_mods" : calc_mods, "spec_qad" : spec_qad}
    args_fc = [atoms, smatrix]
    args_db = [atoms, db_path]
    # Initialize the displacements with phonopy
    fw1 = Firework(PyTask({"func": fw.initialize_phonopy.name,
                           "args": args_init}),
                   name="initialize_phonopy")
    fw2 = Firework(PyTask({"func": fw.calculate_multiple.name,
                           "kwargs": kwargs_cm,
                          "inputs": ["atom_dicts", "workdirs"]}),
                   name="setup_calcs")
    anal_task_list = []
    anal_task_list.append(PyTask({"func": fw.calc_phonopy_force_constants.name,
                                  "args": args_fc,
                                  "inputs": ["calc_atoms"]}))
    anal_task_list.append(PyTask({"func": fw.calc_phonopy_band_structure.name,
                                  "args": [],
                                  "inputs": ["phonon_dict"]}))
    anal_task_list.append(PyTask({"func": fw.calc_phonopy_dos.name,
                                  "args": [[45,45,45]],
                                  "inputs": ["phonon_dict"]}))
    anal_task_list.append(PyTask({"func": fw.add_phonon_to_db.name,
                                  "args": args_db,
                                  "inputs": ["phonon_dict"]}))
    fw3 = Firework(anal_task_list, name="analysis_and_saving")

    workflow = Workflow([fw1, fw2, fw3], {fw1:[fw2], fw2:[fw3]})
    launchpad.add_wf(workflow)
    rapidfire(launchpad, nlaunches=2)
    qlaunch_remote("rapidfire", maxjobs_queue=250, nlaunches=2, remote_host=args.remote_host,
                   remote_user=args.remote_user, remote_password=args.remote_password,
                   remote_config_dir=["/u/tpurcell/git/hilde/examples/fireworks"],
                   remote_setup=False, reserve=True, gss_auth=not args.no_kerberos)
    rapidfire(launchpad, nlaunches=1)

phonon = db.get_phonon(selection=[("supercell_matrix", "=", smatrix),
                                  ("atoms_hash", "=", atoms_hash),
                                  ("calc_hash", "=", calc_hash),
                                  ("has_fc", "=", True)])

bands, labels = get_bands_and_labels(phonon.primitive)
phonon.set_band_structure(bands)
plt = phonon.plot_band_structure(labels=labels)
plt.ylabel('Frequency [THz]')
plt.savefig('phonon_dispersion.pdf')
