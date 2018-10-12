import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from time import time

from hilde.parsers import read_structure, read_output
from hilde.helpers.supercell import find_cubic_cell, make_supercell
from hilde.helpers import cwd, d2k, get_cubicness, clean_matrix
from hilde.phonopy import phono as ph
from hilde.tasks.calculate import compute_forces_socketio, calculate
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator
from ase.dft.kpoints import get_cellinfo, special_paths, bandpath
from ase.io import Trajectory

atoms = read_structure('si.in')
vol = atoms.get_volume()

# aims
command = 'orterun -n 4 aims.x'
species_dir = '/home/knoop/FHIaims/aimsfiles/species_defaults/light'
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
relax_folder = Path().cwd() / f'relax_{atoms.sysname}'
aims = Aims(**{**aims_settings, **relax_settings})

try:
    ratoms = read_output(relax_folder / 'geometry.in.next_step')
except FileNotFoundError:
    calculate(atoms, aims, relax_folder)
    ratoms = read_structure(relax_folder / 'geometry.in.next_step')

n_target = 10
target_size = n_target / len(atoms)
smatrix = find_cubic_cell(cell=ratoms.cell, target_size=target_size)

print(clean_matrix(make_supercell(ratoms, smatrix).cell, eps=1e-8))

phonon, sc, scs = ph.preprocess(ratoms, smatrix.T)

workdir = Path('./{}_{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(
    atoms.sysname, *smatrix.flatten(), vol)).absolute()
workdir.mkdir(exist_ok=True)

k_grid = d2k(sc, 1)
print(f'k_grid supercell: {k_grid}')

sc_settings = {
    'k_grid': k_grid
}
aims = Aims(**{**aims_settings, **sc_settings})

cellinfo = get_cellinfo(ratoms.cell)
path = special_paths[cellinfo.lattice]
bandpath = bandpath(path, atoms.cell)

# reset the phonon object
phonon, sc, scs = ph.preprocess(ratoms, smatrix.T)

print(f'{len(scs)} supercells created.')

socketio_settings = {
    'use_pimd_wrapper': ('localhost', port),
#    'sc_accuracy_forces': 1e-4,
    'k_grid': k_grid
}
aims = Aims(**{**aims_settings, **socketio_settings})
phonopy_log = workdir/'phonopy.log'
traj_file = workdir / 'phonopy.traj'
tmp_dir = workdir / 'socketio'

if traj_file.exists():
    traj = Trajectory(str(traj_file), 'r')
    force_sets = [a.get_forces() for a in traj]
else:
    stime = time()
    force_sets = compute_forces_socketio(cells=scs,
                                         calculator=aims,
                                         port=port,
                                         workdir=tmp_dir,
                                         traj_file=traj_file,
                                         log_file=phonopy_log)
    timing_socket_io = time()-stime
    print(f'.. done in {timing_socket_io:.2f}s')

dos = ph.get_dos(phonon, force_sets=force_sets)
plt.plot(dos[0], dos[1])
plt.savefig(workdir / 'dos.pdf')
