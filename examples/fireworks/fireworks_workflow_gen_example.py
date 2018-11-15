'''Example for generating/submitting a full relaxation->full phonon workflow'''
from argparse import ArgumentParser
from pathlib import Path

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.fireworks_api_adapter.combined_launcher import rapidfire
from hilde.helpers.brillouinzone import get_bands_and_labels, get_sp_points
from hilde.helpers.hash import hash_atoms
from hilde.helpers.utility_functions import get_smatrix
from hilde.parsers import read_structure
from hilde.phonon_db.phonon_db import connect
from hilde.workflows.relax_phonopy import gen_relax_phonopy_wf
from default_aims_settings import (aims_kgrid_conv_settings,
                                   aims_relax_settings_light,
                                   aims_relax_settings_tight,
                                   aims_force_settings)

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
                    help="directory of basis set files for aims on the remote machine",
                    default=None)
parser.add_argument("-rcmd", "--remote_command",
                    help="command used to run aims on the remote machine",
                    default=None)
parser.add_argument("-wd", "--workdir", default=".",
                    help="directory used to calculate the individual atom calculations")
parser.add_argument("--no_kerberos", action="store_true",
                    help="If set do not use gss_api authentication")
parser.add_argument('-g', '--geometry', default='si.in',
                    help="The input geometry file")
args = parser.parse_args()

# Modifications of settings fore remote environment
if args.remote_species_dir:
    aims_kgrid_conv_settings["species_dir"] = args.remote_species_dir
    aims_relax_settings_light["species_dir"] = args.remote_species_dir
    aims_relax_settings_tight["species_dir"] = args.remote_species_dir
    aims_force_settings["species_dir"] = args.remote_species_dir
if args.remote_command:
    aims_kgrid_conv_settings["aims_command"] = args.remote_command
    aims_relax_settings_light["aims_command"] = args.remote_command
    aims_relax_settings_tight["aims_command"] = args.remote_command
    aims_force_settings["aims_command"] = args.remote_command

# Get aims settings
atoms = read_structure(args.geometry)
aims_relax_settings_light['sym_block'] = atoms.symmetry_block
aims_relax_settings_tight['sym_block'] = atoms.symmetry_block
atoms_hash, _ = hash_atoms(atoms)
smatrix = get_smatrix(atoms, n_target=96)
launchpad = LaunchPadHilde()
# Clean up launchpad
launchpad.clean_up_wflow(f"{atoms.get_chemical_formula()}_{atoms_hash}")

wf = gen_relax_phonopy_wf(args.geometry,
                          "postgresql://hilde:hilde@130.183.206.193:5432/phonopy_db",
                          "postgresql://hilde:hilde@localhost:5432/phonopy_db",
                          f"{atoms.get_chemical_formula()}_{atoms_hash}",
                          args.workdir,
                          smatrix,
                          kgrid_conv=aims_kgrid_conv_settings,
                          relax_light=aims_relax_settings_light,
                          relax_tight=aims_relax_settings_tight,
                          force_calc=aims_force_settings,
                          spec_qad_kgrid={"_queueadapter": {"walltime": "00:02:00"}},
                          spec_qad_relax={"_queueadapter": {"nodes": 10, "walltime": "00:30:00"}},
                          spec_qad_forces={"_queueadapter": {"nodes": 32, "walltime": "00:30:00"}})
launchpad.add_wf(wf)
print(args)
# rapidfire(launchpad, launch_dir='.', nlaunches=0, njobs_queue=250, wflow=wf, njobs_block=500,
#           sleep_time=15, reserve=True, remote_host=args.remote_host,
#           remote_user=args.remote_user, remote_password=args.remote_password,
#           remote_config_dir=["/u/tpurcell/git/hilde/examples/fireworks"],
#           gss_auth=not args.no_kerberos)


db_path = "postgresql://hilde:hilde@localhost:5432/phonopy_db"
db = connect(db_path)
phonon = db.get_phonon(1e-5, selection=[("supercell_matrix", "=", smatrix),
                                        ("original_atoms_hash", "=", atoms_hash),
                                        ("has_fc", "=", True),
                                        ("calc_type", "=", "phonons")])
bands, labels = get_bands_and_labels(atoms)
phonon.set_band_structure(bands)
plt = phonon.plot_band_structure(labels=labels)
plt.ylabel('Frequency [THz]')
plt.savefig('phonon_dispersion.pdf')

# Create animation*.ascii files that can be used to animate the phonons.
animation_dir = Path('animation')
animation_dir.mkdir(exist_ok=True)
q_points = get_sp_points(atoms)
for key, value in q_points.items():
    filename = str(animation_dir/"animation_{}.ascii".format(key.lstrip('\\')))
    phonon.write_animation(q_point=value, filename=filename)
print(f'Animation files saved in animation')
