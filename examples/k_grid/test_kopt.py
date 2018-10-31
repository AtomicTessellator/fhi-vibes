from ase.io import read
from ase.calculators.aims import Aims
from hilde.settings import Settings
from hilde.helpers.paths import cwd
from pathlib import Path
from kpointoptimizer import KPointOptimizer
import numpy as np

material = 'aln'
atoms = read('../' + material + '.in', 0, 'aims')

settings = Settings('../../hilde.conf')
species_dir = str(Path(settings.machine.basissetloc) / 'light')
command = settings.machine.aims_command
tmp_dir = Path('./' + material)
tmp_dir.mkdir(parents=True, exist_ok=True)

# Logging
log_settings = {
    'trajectory': str(tmp_dir/'opt.aims.traj'),
    'logfile': tmp_dir/'md.aims.log'}

# DFT
aims_settings = {
    'command': command,
    'species_dir': species_dir,
    'output_level': 'MD_light',
    'relativistic': 'atomic_zora scalar',
    'xc': 'pw-lda',
    'k_grid': [1, 1, 1], # dummy
    'compute_forces': True #, 'compute_analytical_stress': True
}

calc = Aims(**aims_settings)
atoms.calc = calc

def get_forces(atoms):
    return atoms.get_forces()

def get_forces_and_stress(atoms):
    return np.array([*atoms.get_forces().flatten(), *atoms.get_stress().flatten()])

def get_energy(atoms):
    return atoms.get_total_energy() / len(atoms)

def loss_func(arg):
    from numpy import asarray
    return asarray(abs(arg)).max()

opt = KPointOptimizer(atoms,
                      func=get_energy,
                      loss_func=loss_func,
                      dfunc_min=1e-6,
                      trajectory=log_settings['trajectory'],
                      even=True
                      )

with cwd(tmp_dir):
    for _ in opt.irun():
        print(opt.kpts, opt.kpts_density, opt.dfunc)
        pass


