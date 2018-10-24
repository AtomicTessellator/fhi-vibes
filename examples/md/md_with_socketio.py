import sys, os
import numpy as np
from ase.io import read
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import PhononHarmonics
from ase.calculators.calculator import kptdensity2monkhorstpack
from ase import units
# from ase.calculators.lammpsrun import LAMMPS
from hilde.settings import Settings
from hilde.helpers.paths import cwd
from pathlib import Path


# Read settings
port = 27182
settings = Settings('../../hilde.conf')
species_dir = str(Path(settings.machine.basissetloc) / 'light')
command = settings.machine.aims_command
tmp_dir = Path('./tmp')
tmp_dir.mkdir(parents=True, exist_ok=True)

# Read input files
atoms = read('si.conv.in', '0', 'aims')
force_constants = np.loadtxt('force_constants.dat')

# Some parameters
temp = 600*units.kB
k_grid = kptdensity2monkhorstpack(atoms, kptdensity=2)

# Logging
log_settings = {
    'trajectory': str(tmp_dir/'md.aims.traj'),
    'logfile': tmp_dir/'md.aims.log'}

# DFT
aims_settings = {
    'command': command,
    'use_pimd_wrapper': ('localhost', port),
    'species_dir': species_dir,
    'output_level': 'MD_light',
    'xc': 'pw-lda',
    'k_grid': k_grid.tolist(),
    'compute_forces': True
}

calc = Aims(**aims_settings)

# initialize positions + velocities
PhononHarmonics(atoms,
                force_constants,
                temp = temp,
                quantum=False)

# md = Langevin(atoms, temperature=300*units.kB, timestep=.4*units.fs, friction=1e-3, **log_settings)
md = VelocityVerlet(atoms, timestep=4*units.fs, **log_settings)

# run
with SocketIOCalculator(calc, log=sys.stdout, port=port) as calc, cwd(tmp_dir):
    atoms.calc = calc
    md.run(steps=20)
