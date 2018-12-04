# Demonstration of phonon calculations using FireWorks
import os
from ase.build import bulk

from fireworks import Firework, FWorker, LaunchPad, PyTask, FWAction, Workflow

# Minimal hilde inputs to make dictionary conversion easier
from hilde.helpers.hash import hash_atoms_and_calc

import numpy as np
from pathlib import Path
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from hilde.phonon_db.row import phonon_to_dict, PhononRow
from hilde.phonon_db.phonon_db import connect as connect_ph
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import generate_firework

# Combined local/remote queue launching
from hilde.fireworks_api_adapter.combined_launcher import rapidfire
from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
# Use ASE aims calculator
from ase.calculators.aims import Aims

mod_name = __name__


# Aims calculator settings
aims_settings = {
    "xc": "pw-lda",
    "relativistic": "atomic_zora scalar",
    "sc_accuracy_rho": 1E-6,
    "sc_accuracy_forces": 1E-3,
    "sc_iter_limit": 50,
    "mixer": "pulay",
    "n_max_pulay": 10,
    "charge_mix_param": 0.3,
    "k_grid": [4, 4, 4],
    "output_level": "MD_light",
    # Add it but should be changed based on the remote machine's settings
    "species_dir": "~/light",
}

# Use diamond silicon as a test material
si = bulk('Si', 'diamond')
si.calc = Aims(**aims_settings)
si_hash, calc_hash = hash_atoms_and_calc(si)

settings = Settings()
workdir = f"/u/{settings.fireworks.remote_user}/.fireworks/Si/"

print(workdir)
# Initialize phononpy
smatrix = np.array([-1, 1, 1, 1, -1, 1, 1, 1, -1]).reshape(3,3)
db_path = os.getcwd() + "/test_ph.db"

q_spec = {
    # Submission script changes are controled by the _queueadapter dictionary
    "_queueadapter":{
        # Keys are the same that you define in "my_qadapter.yaml"
        "walltime": "00:05:00",
        "nodes": 2
    }
}

fw_settings = {
    "serial": True,
    "fw_name": "phonons",
    "fw_spec": None,
    "out_spec": "relaxed_atoms",
    "mod_spec_add": "calculated_forces",
    "spec": q_spec
}

kwargs_init = {"supercell_matrix": smatrix, "displacement": 0.01}
kwargs_init_fw_out = {"workdir": workdir}

init_fw = generate_firework(
    "hilde.phonopy.workflow.initialize_phonopy_attach_calc",
    "hilde.tasks.fireworks.fw_action_outs.add_phonopy_force_calcs",
    kwargs_init,
    si,
    si.calc,
    func_fw_out_kwargs=kwargs_init_fw_out,
    atoms_calc_from_spec=False,
    fw_settings=fw_settings,
)

kwargs = {
    "fireworks": True,
    "db_path": db_path,
    "original_atom_hash": si_hash,
    "workdir": str(Path("Si").absolute()),
    "displacement": 0.01,
}

anal_fw = generate_firework(
    "hilde.phonopy.postprocess.postprocess",
    "hilde.tasks.fireworks.fw_action_outs.fireworks_no_mods_gen_function",
    args=[],
    inputs=["phonon", fw_settings["mod_spec_add"]],
    func_kwargs=kwargs,
    fw_settings=fw_settings,
)
lp = LaunchPadHilde.from_file(str(Path.home() / ".fireworks/my_launchpad.yaml"))

wf = Workflow([init_fw, anal_fw], {init_fw: [anal_fw]})
lp.add_wf(wf)

rapidfire(lp, FWorker(), strm_lvl="INFO", reserve=True, gss_auth=True)

#Access the database to check the results
db = connect_ph(db_path)
print(si_hash, smatrix)
row = list(db.select(selection=[("original_atom_hash", "=", si_hash),
                                ("sc_matrix_2", "=", smatrix)],
                     columns=["fc_2"]))[0]
print(f"The force constants are:\n{row.get('fc_2')}")