""" An example of how to use FireWorks in conjunction with HilDe"""
import os
from pathlib import Path

import numpy as np

from ase.build import bulk
from ase.calculators.emt import EMT

from hilde.fireworks.launchpad import LaunchPadHilde
from hilde.fireworks.rocket_launcher import rapidfire
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.helpers.paths import cwd
from hilde.fireworks.workflows.workflow_generator import (
    generate_firework,
    get_phonon_analysis_task,
)
from hilde.phonopy.postprocess import postprocess

from fireworks import Workflow

db_name = os.getcwd() + "/test.db"
print(f"database: {db_name}")

atoms = bulk("Ni", "fcc", a=3.5)
atoms.set_calculator(EMT())
calc = atoms.calc
workdir = str(Path("Ni_ex").absolute())

smatrix = np.array([[-2, 2, 2], [2, -2, 2], [2, 2, -2]])
atoms_hash, calc_hash = hash_atoms_and_calc(atoms)

fw_settings = {
    "serial": True,
    "fw_name": "phonons",
    "fw_spec": None,
    "out_spec": "relaxed_atoms",
    "mod_spec_add": "ph_forces",
    "spec": {},
}

# set up the LaunchPad and reset it
kwargs_init = {"supercell_matrix": smatrix, "displacement": 0.03, "serial": True}
kwargs_init_fw_out = {"workdir": workdir, "serial": True}

init_fw = generate_firework(
    func="hilde.tasks.fireworks.phonopy_phono3py_functions.bootstrap_phonon",
    func_fw_out="hilde.tasks.fireworks.fw_out.phonons.post_init_mult_calcs",
    func_kwargs={"ph_settings": kwargs_init},
    atoms=atoms,
    calc=calc,
    args=[1.0],
    func_fw_out_kwargs={"ph_settings": kwargs_init_fw_out},
    atoms_calc_from_spec=False,
    fw_settings=fw_settings,
)

kwargs = {"fireworks": True, "workdir": workdir + "/analysis"}
task_spec_list = get_phonon_analysis_task(
    "hilde.phonopy.postprocess.postprocess", kwargs, "ph_metadata", "ph_forces", "ph_times"
)

anal_fw = generate_firework(task_spec_list, fw_settings=fw_settings)

lp = LaunchPadHilde(strm_lvl="INFO")
lp.reset("", require_password=False)
wf = Workflow([init_fw, anal_fw], {init_fw: [anal_fw]})
lp.add_wf(wf)

with cwd(workdir + "/fireworks", mkdir=True):
    rapidfire(lp, wflow_id=wf.root_fw_ids, strm_lvl="INFO")

phonon = postprocess(workdir + "/analysis/trajectory.yaml")

phonon.set_mesh(3 * [5])
_, _, frequencies, _ = phonon.get_mesh()
print(f"Highest frequency: {frequencies.max():.3f} THz (Target: [8,10] THz)")
assert 8 < frequencies.max() < 10
