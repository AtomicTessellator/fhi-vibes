from ase.io import read, write
from ase.constraints import UnitCellFilter
from ase.optimize import BFGS, FIRE
from ase.optimize.gpmin.gpmin import GPMin
from ase.build import make_supercell
from ase.calculators.socketio import SocketIOCalculator
from vibes.templates.aims import setup_aims

atoms = read("../si.in", 0, "aims")

spos = atoms.get_scaled_positions()
atoms.cell = atoms.cell * 1.2
atoms.set_scaled_positions(spos)

conv_matrix = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
atoms = make_supercell(atoms, conv_matrix)
atoms.rattle()

atoms.write("geometry.in", "aims")

port = 27182
calc = setup_aims(
    custom_settings={
        "use_pimd_wrapper": ("localhost", port),
        "compute_forces": True,
        "compute_analytical_stress": True,
        # "use_symmetric_forces": True
    },
    workdir="tmp",
)

opt_atoms =  UnitCellFilter(atoms)

optimizer = GPMin

# can be used as soon as MR 998 is merged to ASE
# with SocketIOCalculator(calc, log="socketio.log", port=port) as calc:
with SocketIOCalculator(calc, log=open("socketio.log", 'w'), port=port) as calc:
    atoms.set_calculator(calc)
    opt = optimizer(opt_atoms, logfile="relax.log")
    for _ in opt.irun(fmax=0.01, steps=20):
        print(atoms.get_scaled_positions())
        print(atoms.cell)
        print()
        print(opt_atoms.get_forces())
        print(abs(atoms.get_forces()).max())
        print(abs(opt_atoms.get_forces()).max())
        print()

atoms.write("geometry.out", "aims")
