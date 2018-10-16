from glob import glob
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import shutil

from ase.io import read, write

from hilde.helpers.brillouinzone import get_bands, get_sp_points
from hilde.materials_fp.material_fingerprint import *
from hilde.parsers import read_structure, read_output
from hilde.phonopy import phono as ph
from hilde.settings import Settings
from hilde.tasks.calculate import calculate_multiple
from hilde.templates.aims import setup_aims

def make_workdir(smatrix, vol):
    wd = Path('./Si_{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(*smatrix.flatten(), vol))
    wd.mkdir(exist_ok=True)
    return wd
settings = Settings('../../hilde.conf')
aims_settings = {
    'command': settings.machine.aims_command,
    'species_dir': str(Path(settings.machine.basissetloc) / 'light'),
    'output_level': 'MD_light',
    'relativistic': 'atomic_zora scalar',
    'xc': 'pw-lda',
    'k_grid': 3 * [2],
    "sc_accuracy_rho" : 0.0001,
    "sc_accuracy_forces" : 0.0005
}

calc = setup_aims(aims_settings)

# Electronic Modes
kGridDirs = glob("k_grid_conv/*")
kGridDirs.sort(key=lambda s: float(s.split("_")[-1] ) )
kpoints = { "L"      : np.array([0.500, 0.500, 0.500]),
            "\Gamma" : np.array([0.000, 0.000, 0.000]),
            "X"      : np.array([0.000, 0.500, 0.500]),
            "W"      : np.array([0.250, 0.500, 0.750]),
            "K"      : np.array([0.375, 0.375, 0.750])
          }

elecDs = [ DOSFingerprint(True, False, spectra_files=[s + '/KS_DOS_total.dat'], min_e=-15.0, max_e=10.0, nbins=100) for s in kGridDirs ]
elecBs = [ BandStructureFingerprint(True, True, kpoints=kpoints, spectra_files=[s + '/band1001.out', s + '/band1002.out', s + '/band1003.out', s + '/band1004.out'], min_e=-1800.0, max_e=20.0, nbins=20 ) for s in kGridDirs ]

print("DOSFingerprint Convergence")
for ii in range(len(elecDs)):
    print(kGridDirs[ii].split("_")[-1], elecDs[-1].scalar_product(elecDs[ii], 1, 'All'))
print("BandStructureFingerprint Convergence")
for ii in range(len(elecBs)):
    print(kGridDirs[ii].split("_")[-1], elecBs[-1].scalar_product(elecBs[ii], 0, 'All'))

# Phonons
# Set up the phonopy objects
atoms = read_structure('../phonons/geometry.in')
vol = atoms.get_volume()
smatrix = np.array([[-1,  1,  1],
                    [ 1, -1,  1],
                    [ 1,  1, -1]])

# A series of super cell matrices
smatrices = [a * smatrix for a in range(1,3)]
phononCalcs = [ph.preprocess(atoms, smatrix)+(make_workdir(smatrix, vol),) for smatrix in smatrices]

# Calculate the Forces
fp_list = []
for phonon, sc, scs, wd in phononCalcs:
    scs = calculate_multiple(scs, calc, wd, force=True)
    phonon.set_forces([sc.get_forces() for sc in scs])
    phonon.produce_force_constants()

    bands = get_bands(atoms)
    phonon.set_band_structure(bands)

    fp = get_phonon_bs_fingerprint_phononpy(phonon, get_sp_points(atoms), binning=False)
    fp_list.append( fp )

    q_mesh = [45, 45, 45]
    phonon.set_mesh(q_mesh)
    # We generate the DOS by calling .set_total_DOS on the phonopy object
    phonon.set_total_DOS(freq_pitch=.1, tetrahedron_method=True)
print("phonon BS scalar product")
for ff in fp_list:
    print( scalar_product(ff, fp_list[-1], 0, 0, True) )

phonon_dos_fp = to_dict(get_phonon_dos_fingerprint_phononpy(phononCalcs[0][0], binning=True))
print("Phonon DOS")
print(phonon_dos_fp["DOS"])
