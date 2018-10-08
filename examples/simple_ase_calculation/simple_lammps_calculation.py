from ase.io import read
from ase.calculators.lammpsrun import LAMMPS
from hilde.helpers.paths import cwd
from hilde.tasks.calculate import calculate
from ase.calculators.socketio import SocketIOCalculator, PySocketIO
import os

atoms = read('si.in', 0, 'aims')

# LAMMPS context information
lmp_path = os.getenv("LAMMPS_PATH")
potential = os.path.join(lmp_path, "potentials", "Si.tersoff")
files = [potential]
parameters = {"mass": ["* 1.0"],
              "pair_style": "tersoff",
              "pair_coeff": ['* * ' + potential + ' Si']}

lammps = LAMMPS(parameters=parameters,
                files=files)

# option 1:
atoms.calc = lammps
# REM: using .calculate(atoms) is mandatory for lammmps so far
atoms.calc.calculate(atoms)
print(atoms.get_total_energy())
print(atoms.get_forces())
