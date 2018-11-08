""" run an MD for 2 steps with socketio """

import numpy as np
from ase.io import read
from ase.calculators.socketio import SocketIOCalculator
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import PhononHarmonics
from ase import units
from hilde.templates.aims import setup_aims

# Read input files
atoms = read("Si.in.supercell", "0", "aims")
force_constants = np.loadtxt("force_constants_Si.dat")

# use template aims calculator (taking settings from hilde.cfg) and custom settings
port = 1234
custom_settings = {
    "use_pimd_wrapper": ("localhost", port),
    "sc_accuracy_rho": 1e-4,
    "compute_forces": True,
}

calc = setup_aims(custom_settings, workdir="aims_tmp")

# initialize positions + velocities at 100K
PhononHarmonics(atoms, force_constants, temp=100 * units.kB, quantum=False)

md = VelocityVerlet(
    atoms, timestep=4 * units.fs, trajectory="md.aims.traj", logfile="md.aims.log"
)

# run
with SocketIOCalculator(calc, log="socketio.log", port=port) as calc:
    atoms.calc = calc
    md.run(steps=2)
