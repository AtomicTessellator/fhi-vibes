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

n_target = 10
sc, smatrix = make_cubic_supercell(ratoms, target_size=n_target)

print(sc.cell)

phonon, sc, scs = ph.preprocess(ratoms, smatrix.T)

workdir = base_folder / '{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(
    *smatrix.flatten(), vol)
workdir.mkdir(exist_ok=True, parents=True)

sc.write(workdir / 'supercell.in')

k_grid = d2k(sc, 1)
print(f'k_grid supercell: {k_grid}')

# set the phonon object
phonon, sc, scs = ph.preprocess(ratoms, smatrix.T)

print(f'{len(scs)} supercells created.')

socketio_settings = {
    'use_pimd_wrapper': ('localhost', port),
    'k_grid': k_grid
}

aims = Aims(**{**aims_settings, **socketio_settings})

phonopy_log = workdir/'phonopy.log'
trajectory = workdir / 'phonopy.traj'
tmp_dir = workdir / 'socketio'

stime = time()
cells_computed = calculate_multiple_socketio(cells=scs,
                                             calculator=aims,
                                             workdir=tmp_dir,
                                             port=port,
                                             trajectory=trajectory,
                                             force=False,
                                             log_file=phonopy_log)
timing_socket_io = time() - stime
print(f'.. done in {timing_socket_io:.2f}s')

force_sets = [atoms.get_forces() for atoms in cells_computed]

dos = ph.get_dos(phonon, force_sets=force_sets)
plt.plot(dos[0], dos[1])
plt.savefig(workdir / 'dos.pdf')
