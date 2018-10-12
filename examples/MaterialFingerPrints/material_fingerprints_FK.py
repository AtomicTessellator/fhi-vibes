""" Compute the phonon fingerprints for supercells of different size """

from pathlib import Path
import numpy as np
from hilde.parsers import read_structure
from hilde.phonopy import phono as ph
from hilde.tasks import compute_forces
from ase.dft.kpoints import get_cellinfo
from hilde.templates.lammps import setup_lammps_si
from hilde.materials_fp.MaterialsFingerprints import get_phonon_bs_fingerprint_phononpy
from hilde.helpers.supercell import make_cubic_supercell

atoms = read_structure('si.in')

# Space group information
cellinfo = get_cellinfo(atoms.cell)
special_points = cellinfo.special_points

# Calculator setup

# conventional supercell matrix
cmatrix = np.array([[-1,  1,  1],
                    [ 1, -1,  1],
                    [ 1,  1, -1]])

# run phonon calculation for several supercell sizes and compute fingerprints
fps = []
n_atoms = []
for nn in [8, 64, 128, 216]:

    # smatrix = a * cmatrix
    supercell, smatrix = make_cubic_supercell(atoms, nn)

    n_a = len(supercell)
    print(f'compute for {n_a} atoms')
    n_atoms.append(n_a)

    phonon, sc, scs = ph.preprocess(atoms, smatrix.T)

    workdir = Path('./Si_{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(
        *smatrix.flatten(), atoms.get_volume())).absolute()
    workdir.mkdir(exist_ok=True)

    lammps = setup_lammps_si(workdir)

    force_sets = compute_forces(scs, lammps, workdir)
    phonon.produce_force_constants(force_sets)

    fp = get_phonon_bs_fingerprint_phononpy(phonon, special_points)
    fps.append([p[:, 0] for p in fp.values()])

fps = np.asarray(fps)

# Compute difference to largest supercell and choose largest deviation at each
# q point
fp_diffs = abs(fps - fps[-1]).max(axis=2)

print('n_atoms  '  + ' '.join([f'{k:9s}' for k in special_points.keys()]))
for nn, fp in zip(n_atoms, fp_diffs):
    print(f'{nn:4d}: ' + ' '.join([f"{f:9.3e}" for f in fp]))
