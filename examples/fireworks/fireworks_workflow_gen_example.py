from argparse import ArgumentParser
import numpy as np
import os

from fireworks import Firework, LaunchPad, PyTask, Workflow
from fireworks.core.rocket_launcher import rapidfire, launch_rocket

from hilde.fireworks_api_adapter.qlaunch_remote import qlaunch_remote
from hilde.helpers.brillouinzone import get_bands_and_labels
from hilde.helpers.paths import cwd
from hilde.helpers.utility_functions import get_smatrix, setup_workdir
from hilde.parsers import read_structure
from hilde.phonon_db.phonon_db import connect
from hilde.workflows.relax_phonopy import gen_relax_phonopy_wf

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
parser.add_argument("--no_kerberos", action="store_true",
                    help="If set do not use gss_api authentication")
args = parser.parse_args()
atoms = read_structure('si.in')
smatrix = get_smatrix(atoms, n_target=64)
launchpad = LaunchPad()
try:
    query = {"name": "example_SI_wf_fill", "state": "COMPLETED"}
    wf_ids = launchpad.get_wf_ids(query=query, limit=100)
    for wf_id in wf_ids:
        launchpad.delete_wf(wf_id)
except:
  pass
try:
    query = {"name": "example_SI_wf_fill", "state": "FIZZLED"}
    wf_ids = launchpad.get_wf_ids(query=query, limit=100)
    for wf_id in wf_ids:
        launchpad.delete_wf(wf_id)
except:
    pass
try:
    query = {"name": "example_SI_wf_fill", "state": "READY"}
    wf_ids = launchpad.get_wf_ids(query=query, limit=100)
    for wf_id in wf_ids:
        launchpad.delete_wf(wf_id)
except:
    pass
try:
    query = {"name": "example_SI_wf_fill", "state": "RESERVED"}
    wf_ids = launchpad.get_wf_ids(query=query, limit=100)
    for wf_id in wf_ids:
        launchpad.delete_wf(wf_id)
except:
    pass
wf = gen_relax_phonopy_wf("si.in",
                          "postgresql://hilde:hilde@130.183.206.193:5432/phonopy_db",
                          "example_SI_wf_fill",
                          "/u/tpurcell/git/hilde/examples/fireworks/test_wf",
                          "atoms_cur",
                          smatrix,
                          symprec=1e-2,
                          spec_qad_kgrid={"'_queueadapter'" : {"walltime" : "00:10:00"} },
                          spec_qad_relax={"'_queueadapter'" : {"walltime" : "00:10:00"} },
                          spec_qad_forces={"'_queueadapter'" : {"nodes" : 3} })

launchpad.add_wf(wf)
qlaunch_remote("rapidfire", maxjobs_queue=250, nlaunches=0, remote_host=args.remote_host,
               remote_user=args.remote_user, remote_password=args.remote_password,
               remote_config_dir=["/u/tpurcell/git/hilde/examples/fireworks"], reserve=True,
               gss_auth=not args.no_kerberos)

db_path = "postgresql://hilde:hilde@localhost:5432/phonopy_db"
db = connect(db_path)
phonon = db.get_phonon(1e-2, selection=[("supercell_matrix", "=", smatrix),
                                  ("numbers", "=", atoms.numbers),
                                  ("has_fc", "=", True)])
bands, labels = get_bands_and_labels(atoms)
phonon.set_band_structure(bands)
plt = phonon.plot_band_structure(labels=labels)
plt.ylabel('Frequency [THz]')
plt.savefig('phonon_dispersion.pdf')
