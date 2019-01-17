import numpy as np

from ase.io import read, Trajectory
from ase.constraints import UnitCellFilter

# from ase.optimize import BFGSLineSearch as BFGS
from ase.optimize import BFGS as BFGS
from ase.calculators.socketio import SocketIOCalculator
from ase.calculators.aims import Aims

# read input geometry that should be relaxed
atoms = read("../si.in", 0, "aims")

# choose a name for an aims working directory
workdir = "aims"
port = 12345

# user settings (adjust!)
usr_settings = {
    "aims_command": "orterun aims.x",
    "species_dir": "/home/knoop/FHIaims/aimsfiles/species_defaults/light",
    "output_level": "MD_light",
    "outfilename": "aims.out",
}

# dft settings (adjust!)
dft_settings = {
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "vdw_correction_hirshfeld": True,
    "sc_accuracy_rho": 1e-6,
    "k_grid": [4, 4, 4],
    "compute_forces": True,
    "compute_analytical_stress": True,
    "use_symmetric_forces": True,
}

# BFGSLinesearch settings
bfgs_settings = {"logfile": "relax.log", "trajectory": "relax.traj"}

# some more settings that can stay untouched
aux_settings = {"label": workdir, "use_pimd_wrapper": ("localhost", port)}

# create Aims calculator from the settings
calc = Aims(**usr_settings, **dft_settings, **aux_settings)

# Create a filter to expose the lattice degrees of freedom to the optimizer
opt_atoms = UnitCellFilter(atoms)

# choose an optimizer
optimizer = BFGS


# run the optimization and write out every step
with SocketIOCalculator(calc, log="socketio.log", port=port) as calc:
    atoms.set_calculator(calc)
    # create the optimizer for the filtered atoms object
    opt = optimizer(opt_atoms, **bfgs_settings)
    for ii, _ in enumerate(opt.irun(fmax=0.01, steps=20)):
        # apply external pressure
        stress = atoms.get_stress() + [0.1, 0.1, 0.1, 0, 0, 0]
        atoms.calc.results['stress'] = stress

# write final result
atoms.write("geometry.out", "aims", scaled=True)

# read the relaxation trajectory
with Trajectory(bfgs_settings["trajectory"], "r") as traj:
    trajectory = [at for at in traj]

# compute deformation:
# new_lattice = ( 1 + deformation ) * old_lattice
deformation = np.linalg.inv(trajectory[0].cell) @ trajectory[-1].cell - np.eye(3)

print(f"Optimization resulted in a deformation of")
print(np.around(deformation, decimals=4))
