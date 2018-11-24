''' An example of how to use FireWorks in conjunction with HilDe'''
import importlib as il
import os
from ase.build import bulk
from ase.calculators.emt import EMT

from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.fireworks_api_adapter.rocket_launcher import rapidfire
from hilde.helpers.hash import hash_atoms
from hilde.helpers.paths import cwd
from hilde.helpers.utility_functions import get_smatrix
from hilde.phonon_db.phonon_db import connect
from hilde.structure.structure import pAtoms
from hilde.workflows.gen_phono3py_fw import gen_initialize_phono3py_fw, gen_analyze_phono3py_fw

from fireworks import Workflow
ph3 = il.import_module('hilde.phono3py.phono3')

db_name = (os.getcwd() + '/test.json')
print(f'database: {db_name}')

atoms = pAtoms(bulk('Al'))
atoms.set_calculator(EMT())
calc = atoms.calc
workdir = "Al"

smatrix = get_smatrix(atoms, n_target=32)
q_mesh = [5, 5, 5]

phono3py_settings = {
    'fc2_supercell_matrix': smatrix,
    'fc3_supercell_matrix': smatrix,
    'cutoff_pair_distance': 3,
    'log_level': 0
}

atoms_hash, calc_hash = hash_atoms(atoms)

db = connect(db_name)

# set up the LaunchPad and reset it
launchpad = LaunchPadHilde()
launchpad.clean_up_wflow("Ex_WF_Ni")
fw_get_forces = gen_initialize_phono3py_fw(atoms,
                                           phono3py_settings,
                                           workdir,
                                           calc=calc,
                                           symprec=1e-5,
                                           name="ex_Al_get_force")
fw_analyze = gen_analyze_phono3py_fw(atoms,
                                     db_name,
                                     phono3py_settings)

workflow = Workflow([fw_get_forces, fw_analyze], {fw_get_forces: fw_analyze},
                    name="Al_example")
launchpad.add_wf(workflow)

with cwd(workdir + '/fireworks', mkdir=True):
    rapidfire(launchpad, wflow_id=workflow.root_fw_ids)

phonon3 = db.get_phonon3(selection=[("sc_matrix_3", "=", smatrix),
                                    ("sc_matrix_2", "=", smatrix),
                                    ("atoms_hash", "=", atoms_hash),
                                    ("calc_hash", "=", calc_hash),
                                    ("has_fc3", "=", True),
                                    ("has_fc2", "=", True),
                                    ("calc_type", "=", "phonons")])
fc2 = phonon3.get_fc2()
fc3 = phonon3.get_fc3()
phonon3 = ph3.prepare_phono3py(atoms,
                               **phono3py_settings,
                               fc2=fc2,
                               fc3=fc3,
                               q_mesh=q_mesh)
phonon3.run_thermal_conductivity(write_kappa=True, temperatures=[300])
print(f"Thermal conductivity of the material is: {phonon3.get_thermal_conductivity().get_kappa()[0][0][0]}")
assert 9 < phonon3.get_thermal_conductivity().get_kappa()[0][0][0] < 10

