import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from time import time

from hilde.parsers import read_structure, read_output
from hilde.helpers.supercell import make_cubic_supercell
from hilde.helpers import cwd, d2k, get_cubicness, clean_matrix
from hilde.phonopy import phono as ph
from hilde.tasks.calculate import calculate, calculate_multiple_socketio
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator
from ase.dft.kpoints import get_cellinfo, special_paths, bandpath
from hilde.settings import Settings

atoms = read_structure('../gan.in')
vol = atoms.get_volume()

base_folder = Path(f'{atoms.sysname}').absolute()

# aims
settings = Settings('../../hilde.conf')

command = settings.machine.aims_command
species_dir = str(Path(settings.machine.basissetloc) / 'light')
port = 31415

k_grid = d2k(atoms, 1)
print(f'k_grid primitive cell: {k_grid}')

aims_settings = {
    'command': command,
    'species_dir': species_dir,
    'output_level': 'MD_light',
    'relativistic': 'atomic_zora scalar',
    'xc': 'pw-lda',
    'k_grid': k_grid,
    'sc_accuracy_rho': 1e-6,
    'compute_forces': True
}

relax_settings = {
    'relax_geometry': 'lattice_trm 1e-3',
    'relax_unit_cell': 'full',
    'sc_accuracy_forces': 1e-4,
    'use_symmetric_forces': True
}

relax_folder = base_folder / f'relax_{vol:.3f}'

aims = Aims(**{**aims_settings, **relax_settings})

try:
    ratoms = read_structure(relax_folder / 'geometry.in.next_step')
except FileNotFoundError:
    calculate(atoms, aims, relax_folder)
    ratoms = read_structure(relax_folder / 'geometry.in.next_step')
