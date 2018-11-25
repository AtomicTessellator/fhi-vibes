''' An example of how to use FireWorks in conjunction with HilDe'''
import os
from ase.build import bulk
from ase.calculators.emt import EMT
from pathlib import Path

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.fireworks_api_adapter.rocket_launcher import rapidfire
from hilde.helpers.hash import hash_atoms
from hilde.helpers.paths import cwd
from hilde.helpers.utility_functions import get_smatrix
from hilde.phonon_db.phonon_db import connect
from hilde.structure.structure import pAtoms
from hilde.tasks.fireworks.general_py_task import generate_firework
from hilde.workflows.gen_phonopy_fw import (
    gen_initialize_phonopy_fw,
    gen_analyze_phonopy_fw,
)
from fireworks import Workflow

db_name = (os.getcwd() + '/test.db')
print(f'database: {db_name}')

atoms = pAtoms(bulk('Ni', 'fcc', a=3.5))
atoms.set_calculator(EMT())
calc = atoms.calc
workdir = str(Path("Ni_ex").absolute())

smatrix = get_smatrix(atoms, n_target=32)
atoms_hash, calc_hash = hash_atoms(atoms)

fw_settings = {
    "serial": True,
    "fw_name": "phonons",
    "fw_spec": None,
    "out_spec": "relaxed_atoms",
    "mod_spec_add": "calculated_forces",
    "spec": {}
}

# set up the LaunchPad and reset it

kwargs_init = {"supercell_matrix": smatrix, "displacement": 0.01}
kwargs_init_fw_out = {"workdir": workdir}

init_fw = generate_firework(
    "hilde.phonopy.workflow.initialize_phonopy_attach_calc",
    "hilde.tasks.fireworks.fw_action_outs.add_phonopy_force_calcs",
    kwargs_init,
    atoms,
    calc,
    func_fw_out_kwargs=kwargs_init_fw_out,
    atoms_calc_from_spec=False,
    fw_settings=fw_settings,
)

kwargs = {
    "fireworks": True,
    "db_path": db_name,
    "original_atom_hash": atoms_hash,
    "workdir": workdir,
    "displacement": 0.01,
    "atoms_hash": atoms_hash,
    "calc_hash": calc_hash,
    "has_fc2": True,
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

with cwd(workdir + '/fireworks', mkdir=True):
    rapidfire(lp, wflow_id=wf.root_fw_ids)

db = connect(db_name)

phonon = db.get_phonon(selection=[("sc_matrix_2", "=", smatrix),
                                  ("atoms_hash", "=", atoms_hash),
                                  ("calc_hash", "=", calc_hash),
                                  ("has_fc2", "=", True),
                                  ("calc_type", "=", "phonons")])

phonon.set_mesh(3 * [3])
_, _, frequencies, _ = phonon.get_mesh()
print(f'Highest frequency: {frequencies.max():.3f} THz (Target: [8,10] THz)')
assert 8 < frequencies.max() < 10
