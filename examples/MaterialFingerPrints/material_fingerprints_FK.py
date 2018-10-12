from pathlib import Path
import numpy as np
from hilde.parsers import read_structure
from hilde.phonopy import phono as ph
from hilde.tasks import compute_forces
from hilde.helpers.k_grid import d2k
from ase.dft.kpoints import get_cellinfo
from hilde.templates.lammps import setup_lammps_si
from hilde.materials_fp.MaterialsFingerprints import get_phonon_bs_fingerprint_phononpy

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
for a in range(1, 2):
    smatrix = cmatrix * a

    phonon, sc, scs = ph.preprocess(atoms, smatrix.T)

    workdir = Path('./Si_{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(
        *smatrix.flatten(), atoms.get_volume())).absolute()
    workdir.mkdir(exist_ok=True)

    lammps = setup_lammps_si(workdir)

    force_sets = compute_forces(scs, lammps, workdir)
    phonon.produce_force_constants(force_sets)

    fp = get_phonon_bs_fingerprint_phononpy(phonon, special_points)
    fps.append(fp)

print(fps)
