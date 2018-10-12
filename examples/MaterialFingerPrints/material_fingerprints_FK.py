from pathlib import Path
import numpy as np
from hilde.parsers import read_structure
from hilde.phonopy import phono as ph
from hilde.tasks import compute_forces
from hilde.helpers.k_grid import d2k
from ase.dft.kpoints import get_cellinfo
from ase.calculators.aims import Aims
from hilde.materials_fp.MaterialsFingerprints import get_phonon_bs_fingerprint_phononpy

atoms = read_structure('si.in')

# Space group information
cellinfo = get_cellinfo(atoms.cell)
special_points = cellinfo.special_points

# Calculator setup
command = 'orterun -n 4 aims.x'
species_dir = '/home/knoop/FHIaims/aimsfiles/species_defaults/light'
aims_settings = {
    'command': command,
    'species_dir': species_dir,
    'output_level': 'MD_light',
    'relativistic': 'atomic_zora scalar',
    'xc': 'pw-lda',
    'k_grid': [2, 2, 2],
    'sc_accuracy_rho': 1e-5,
    'compute_forces': True}

# conventional supercell matrix
cmatrix = np.array([[-1,  1,  1],
                    [ 1, -1,  1],
                    [ 1,  1, -1]])

# run phonon calculation for several supercell sizes and compute fingerprints
fps = []
for a in range(1, 2):
    smatrix = cmatrix * a

    phonon, sc, scs = ph.preprocess(atoms, smatrix.T)

    # Adjust k_grid according to supercell size
    sc_settings = {'k_grid': d2k(sc, 1)}
    aims = Aims(**{**aims_settings, **sc_settings})

    workdir = Path('./Si_{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(
        *smatrix.flatten(), atoms.get_volume())).absolute()
    workdir.mkdir(exist_ok=True)

    force_sets = compute_forces(scs, aims, workdir)
    phonon.produce_force_constants(force_sets)

    fp = get_phonon_bs_fingerprint_phononpy(phonon, special_points)
    fps.append(fp)

print(fps)