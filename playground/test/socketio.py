import sys, os
from ase.build import molecule
from ase.calculators.aims import Aims
from ase.optimize import BFGS
from ase.calculators.socketio import SocketIOCalculator
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.emt import EMT
from ase.io import read
from pathlib import Path
from playground.settings import Settings
from playground.helpers.paths import cd

port = 27182
st = Settings('hilde.conf')
species_dir = str(Path(st.machine.basissetloc) / 'light')
command = st.machine.aims_command
tmp_dir = Path('./test/tmp')
tmp_dir.mkdir(parents=True, exist_ok=True)

# Si parameters
# parameters = {"mass": ["* 1.0"],
#               "pair_style": "tersoff",
#               "pair_coeff": ['* * ' + potential + ' Si']}

# H2O parameters:

# calc = EMT()

atoms = read('./test/h2o.xyz', '0', 'xyz')
atoms.rattle(stdev=0.1)

# MaxwellBoltzmannDistribution(atoms, 0.5 * 300 * units.kB, force_temp=True)

calc = Aims(command=command,
            use_pimd_wrapper=('localhost', port),
            compute_forces=True,
            xc='LDA',
            species_dir=species_dir,
            output_level='MD_light')

md = Langevin(atoms, temperature=300*units.kB, timestep=.5*units.fs, friction=1e-3,
              trajectory=str(tmp_dir/'md.aims.traj'), logfile=tmp_dir/'md.aims.log')

with SocketIOCalculator(calc, log=sys.stdout, port=port) as calc, cd(tmp_dir):
    atoms.calc = calc
    md.run(steps=10)
