''' An example of how to use FireWorks in conjunction with HilDe'''
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
from hilde.workflows.relax_phonopy import gen_initialize_phonopy_fw, gen_analyze_phonopy_fw

from fireworks import Workflow

db_name = (os.getcwd() + '/test.json')
print(f'database: {db_name}')

atoms = pAtoms(bulk('Ni', 'fcc', a=3.5))
atoms.set_calculator(EMT())
calc = atoms.calc
workdir = "Ni_ex"

smatrix = get_smatrix(atoms, n_target=32)
atoms_hash, calc_hash = hash_atoms(atoms)

db = connect(db_name)
has_fc = False
found = False

# set up the LaunchPad and reset it
launchpad = LaunchPadHilde()
launchpad.clean_up_wflow("Ex_WF_Ni")
fw_get_forces = gen_initialize_phonopy_fw(atoms,
                                          smatrix,
                                          workdir,
                                          calc=calc,
                                          symprec=1e-5,
                                          name="ex_Ni_get_force")
fw_analyze = gen_analyze_phonopy_fw(atoms,
                                    db_name,
                                    smatrix)

workflow = Workflow([fw_get_forces, fw_analyze], {fw_get_forces: fw_analyze},
                    name="Ex_WF_Ni")
launchpad.add_wf(workflow)

with cwd(workdir + '/fireworks', mkdir=True):
    rapidfire(launchpad, wflow_id=workflow.root_fw_ids)

phonon = db.get_phonon(1e-5, selection=[("supercell_matrix", "=", smatrix),
                                        ("atoms_hash", "=", atoms_hash),
                                        ("calc_hash", "=", calc_hash),
                                        ("has_fc", "=", True),
                                        ("calc_type", "=", "phonons")])

phonon.set_mesh(3 * [3])
_, _, frequencies, _ = phonon.get_mesh()
print(f'Highest frequency: {frequencies.max():.3f} THz (Target: [8,10] THz)')
assert 8 < frequencies.max() < 10
