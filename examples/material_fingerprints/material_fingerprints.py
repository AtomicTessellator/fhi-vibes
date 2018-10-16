'''Example use of material fingerprints'''
from glob import glob
from pathlib import Path
import numpy as np

from hilde.helpers.brillouinzone import get_bands, get_sp_points
from hilde.materials_fp.material_fingerprint import get_phonon_bs_fingerprint_phononpy
from hilde.materials_fp.material_fingerprint import get_phonon_dos_fingerprint_phononpy
from hilde.materials_fp.material_fingerprint import scalar_product
from hilde.materials_fp.material_fingerprint import to_dict
from hilde.materials_fp.material_fingerprint import DOSFingerprint
from hilde.materials_fp.material_fingerprint import BandStructureFingerprint
from hilde.parsers import read_structure
from hilde.phonopy import phono as ph
from hilde.settings import Settings
from hilde.tasks.calculate import calculate_multiple
from hilde.templates.aims import setup_aims

def make_workdir(smat, volume):
    '''Function that creates a working directory from the smat and cell volume'''
    workdirec = Path('./Si_{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(*smat.flatten(), volume))
    workdirec.mkdir(exist_ok=True)
    return workdirec

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
k_grid_dirs = glob("k_grid_conv/*")
k_grid_dirs.sort(key=lambda s: float(s.split("_")[-1]))
kpoints = {"l": np.array([0.500, 0.500, 0.500]),
           "G": np.array([0.000, 0.000, 0.000]),
           "x": np.array([0.000, 0.500, 0.500]),
           "w": np.array([0.250, 0.500, 0.750]),
           "k": np.array([0.375, 0.375, 0.750])
          }

elec_d = [DOSFingerprint(True,
                         False,
                         spectra_files=[s + '/KS_DOS_total.dat'],
                         min_e=-15.0,
                         max_e=10.0,
                         nbins=100
                        )
          for s in k_grid_dirs
         ]

elec_b = [BandStructureFingerprint(True,
                                   True,
                                   kpoints=kpoints,
                                   spectra_files=[s + '/band1001.out',
                                                  s + '/band1002.out',
                                                  s + '/band1003.out',
                                                  s + '/band1004.out'
                                                 ],
                                   min_e=-1800.0,
                                   max_e=20.0,
                                   nbins=20
                                  )
          for s in k_grid_dirs
         ]

print("DOSFingerprint Convergence")
for ii, ed in enumerate(elec_d):
    print(k_grid_dirs[ii].split("_")[-1], elec_d[-1].scalar_product(ed, 1, 'All'))
print("BandStructureFingerprint Convergence")
for ii, eb in enumerate(elec_b):
    print(k_grid_dirs[ii].split("_")[-1], elec_b[-1].scalar_product(eb, 0, 'All'))

# Phonons
# Set up the phonopy objects
atoms = read_structure('../phonons/geometry.in')
vol = atoms.get_volume()
smatrix = np.array([[-1, 1, 1],
                    [1, -1, 1],
                    [1, 1, -1]])

# A series of super cell matrices
smatrices = [a * smatrix for a in range(1, 3)]
phonon_calcs = [ph.preprocess(atoms, sm)+(make_workdir(sm, vol),) for sm in smatrices]

# Calculate the Forces
fp_list = []
for phonon, sc, scs, wd in phonon_calcs:
    scs = calculate_multiple(scs, calc, wd, force=True)
    phonon.set_forces([sc.get_forces() for sc in scs])
    phonon.produce_force_constants()

    bands = get_bands(atoms)
    phonon.set_band_structure(bands)

    fp = get_phonon_bs_fingerprint_phononpy(phonon, get_sp_points(atoms), binning=False)
    fp_list.append(fp)

    q_mesh = [45, 45, 45]
    phonon.set_mesh(q_mesh)
    # We generate the DOS by calling .set_total_DOS on the phonopy object
    phonon.set_total_DOS(freq_pitch=.1, tetrahedron_method=True)
print("phonon BS scalar product")
for ff in fp_list:
    print(scalar_product(ff, fp_list[-1], 0, 0, True))

phonon_dos_fp = to_dict(get_phonon_dos_fingerprint_phononpy(phonon_calcs[0][0], binning=True))
print("Phonon DOS")
print(phonon_dos_fp["DOS"])
