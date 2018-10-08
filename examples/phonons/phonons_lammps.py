import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from time import time

from hilde.parsers import read_aims, read_aims_output
from hilde.helpers.supercell import find_cubic_cell, make_supercell
from hilde.helpers import cwd, d2k, get_cubicness, clean_matrix
from hilde.phonopy import phono as ph
from hilde.tasks.calculate import compute_forces, calculate
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.socketio import SocketIOCalculator
from ase.dft.kpoints import get_cellinfo, special_paths, bandpath
from ase.io import Trajectory

atoms = read_aims('si.in')
vol = atoms.get_volume()


n_target = 10
target_size = n_target / len(atoms)
smatrix = find_cubic_cell(cell=atoms.cell, target_size=target_size)

print(clean_matrix(make_supercell(atoms, smatrix).cell, eps=1e-8))

workdir = Path('./{}_{}{}{}_{}{}{}_{}{}{}_{:.3f}_lammps'.format(
    atoms.sysname, *smatrix.flatten(), vol)).absolute()
workdir.mkdir(exist_ok=True)

# LAMMPS context information
lmp_path = Path(os.getenv("LAMMPS_PATH"))
potential = str(lmp_path / "potentials" / "Si.tersoff")
files = [potential]
parameters = {"mass": ["* 1.0"],
              "pair_style": "tersoff",
              "pair_coeff": ['* * ' + potential + ' Si']}

# Logging
lammps = LAMMPS(parameters=parameters,
                files=files,
                tmp_dir=workdir / 'lammps')

# set the phonon object
phonon, sc, scs = ph.preprocess(atoms, smatrix.T)
print(f'{len(scs)} supercells created.')

phonopy_log = workdir/'phonopy.log'
traj_file = workdir / 'phonopy.traj'
tmp_dir = workdir

if traj_file.exists():
    traj = Trajectory(str(traj_file), 'r')
    force_sets = [a.get_forces() for a in traj]
else:
    stime = time()
    force_sets = compute_forces(cells=scs,
                                calculator=lammps,
                                workdir=tmp_dir)
    timing_socket_io = time()-stime
    print(f'.. done in {timing_socket_io:.2f}s')

dos = ph.get_dos(phonon, force_sets=force_sets)
plt.plot(dos[0], dos[1])
plt.savefig(workdir / f'dos.pdf')

cellinfo = get_cellinfo(atoms.cell)
path = special_paths[cellinfo.lattice]
bands = bandpath(path, atoms.cell)

plt.figure()
qp, dd, fr, ev = ph.get_bandstructure(phonon, bands)
_ = plt.plot(dd[0], fr[0])
plt.savefig(workdir / 'bandstructure.pdf')
exit()

phonon.set_band_structure(list(bands))
print(path)

plt = phonon.plot_band_structure(labels=path.split())
plt.ylabel('Frequency [THz]')
# We save the plot in the working directory.
plt.savefig(str(workdir / 'bandstructure.pdf'))
