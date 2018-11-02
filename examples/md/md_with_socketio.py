""" run an MD with socketio """

import sys
from pathlib import Path
import numpy as np
from ase.io import read
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import PhononHarmonics
from ase import units
from hilde.settings import Settings
from hilde.helpers.paths import cwd


# Read settings
port = 27182
settings = Settings("../../hilde.cfg")
species_dir = str(Path(settings.machine.basissetloc) / "light")
command = settings.machine.aims_command
tmp_dir = Path("./tmp")
tmp_dir.mkdir(parents=True, exist_ok=True)

# Read input files
atoms = read("Si.in.supercell", "0", "aims")
force_constants = np.loadtxt("force_constants_Si.dat")

# Some parameters
temp = 100 * units.kB

# Logging
log_settings = {"trajectory": str(tmp_dir / "md.aims.traj"), "logfile": "md.aims.log"}

# DFT
aims_settings = {
    "command": command,
    "use_pimd_wrapper": ("localhost", port),
    "species_dir": species_dir,
    "output_level": "MD_light",
    "xc": "pw-lda",
    "k_grid": [2, 2, 2],
    "sc_accuracy_rho": 1e-4,
    "compute_forces": True,
}

calc = Aims(**aims_settings)

# initialize positions + velocities
PhononHarmonics(atoms, force_constants, temp=temp, quantum=False)

md = VelocityVerlet(atoms, timestep=4 * units.fs, **log_settings)

# run
with SocketIOCalculator(calc, log=sys.stdout, port=port) as calc, cwd(tmp_dir):
    atoms.calc = calc
    md.run(steps=2)
