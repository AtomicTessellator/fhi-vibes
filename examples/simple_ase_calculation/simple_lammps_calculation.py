from ase.io import read
from ase.calculators.lammpsrun import LAMMPS
from hilde.helpers.paths import cwd
from hilde.tasks.calculate import calculate
from ase.calculators.socketio import SocketIOCalculator, PySocketIO
import os


# Lammps generic
lmp_path = os.getenv("LAMMPS_PATH")

# Silicon
# this is what hilde.templates.lammps.setup_si does
def setup_si():
    atoms = read('si.in', 0, 'aims')
    potential = os.path.join(lmp_path, "potentials", "Si.tersoff")
    files = [potential]
    parameters = {"mass": ["* 1.0"],
                "pair_style": "tersoff",
                "pair_coeff": ['* * ' +  potential + ' Si']}

    return atoms, parameters, files

# GaN
# this is what hilde.templates.lammps.setup_gan does
def setup_gan():
    atoms = read('gan.in', 0, 'aims')
    potential = os.path.join(lmp_path, "potentials", "GaN.tersoff")
    files = [potential]
    parameters = {"mass": ["* 1.0"],
                "pair_style": "tersoff",
                "pair_coeff": ['* * ' +  potential + ' Ga N']}

def main():
    atoms, parameters, files = setup_si()

    lammps = LAMMPS(parameters=parameters,
                    files=files,
                    tmp_dir='./lammps')


    # option 1:
    atoms.calc = lammps
    # REM: using .calculate(atoms) is mandatory for lammmps so far
    atoms.calc.calculate(atoms)
    print(atoms.get_total_energy())
    print(atoms.get_forces())


main()
