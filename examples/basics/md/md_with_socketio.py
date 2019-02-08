""" run an MD for 2 steps with socketio """

from ase.io import read
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

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

# MD settings (adjust!)
md_settings = {"temperature": 600, "timestep": 4}  # K  # fs


# logging settings
log_settings = {"trajectory": "md.traj", "logfile": "md.log"}

# create Aims calculator with socket communication enabled
port = 12345
mandatory_settings = {"use_pimd_wrapper": ("localhost", port), "compute_forces": True}

calc = Aims(**{**mandatory_settings, **dft_settings, **usr_settings})

# Create ASE-MD object, example: VelocityVerlet
md = VelocityVerlet(atoms, timestep=md_settings["timestep"] * units.fs, **log_settings)

# initialize velocities at 600K
MaxwellBoltzmannDistribution(atoms, temp=md_settings["temperature"] * units.kB)


# run
with SocketIOCalculator(calc, log="socketio.log", port=port) as iocalc:
    atoms.calc = iocalc
    md.run(steps=2)

# read trajectory:
from ase.io import Trajectory

with Trajectory(log_settings["trajectory"], "r") as reader:
    trajectory = [atoms for atoms in reader]

# do sth with the atoms objects, e.g.:
# temperatures = [atoms.get_temperature() for atoms in trajectory]
# ...
