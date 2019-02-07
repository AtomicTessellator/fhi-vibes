""" run an MD for 2 steps with socketio """

from ase.io import read
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from ase.io import Trajectory

# Read input files
atoms = read("geometry.in")

# user settings (adjust!)
usr_settings = {
    "aims_command": "orterun aims.x",
    "species_dir": "/home/knoop/FHIaims/aimsfiles/species_defaults/light",
    "output_level": "MD_light",
    "outfilename": "aims.out",
    "label": "calculations",
}


# dft settings (adjust!)
dft_settings = {
    "xc": "pbesol",
    "relativistic": "atomic_zora scalar",
    "sc_accuracy_rho": 1e-6,
    "k_grid": [2, 2, 2],
    "compute_forces": True,
}

# logging settings
log_settings = {"trajectory": "md.traj", "logfile": "md.log"}

# create Aims calculator with socket communication enabled
port = 12345
mandatory_settings = {"use_pimd_wrapper": ("localhost", port), "compute_forces": True}

calc = Aims(**{**mandatory_settings, **dft_settings, **usr_settings})

# initialize velocities at 600K
MaxwellBoltzmannDistribution(atoms, temp=600 * units.kB)

md = VelocityVerlet(atoms, timestep=4 * units.fs, **log_settings)

# run
with SocketIOCalculator(calc, log="socketio.log", port=port) as calc:
    atoms.calc = calc
    md.run(steps=2)

# read trajectory:
with Trajectory(log_settings["trajectory"], "r") as reader:
    trajectory = [atoms for atoms in reader]

# do sth with the atoms objects, e.g.:
# temperatures = [atoms.get_temperature() for atoms in trajectory]
# ...
