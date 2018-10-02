from pathlib import Path
import numpy as np
import shutil
from hilde.parsers import read_aims, read_aims_output
from hilde.phonopy import phono as ph
from hilde.structure.convert import phonopy_to_ASE_atoms
from ase.dft.kpoints import get_cellinfo, special_paths, bandpath
from ase.io import read, write
import matplotlib.pyplot as plt

import sys
from hilde.materials_fp.MaterialsFingerprints import *
import numpy as np
from glob import glob

def get_band(q_start, q_end, n_q=51):
    """ Create a list of n_q equidistant qpoints connecting q_start and q_end """
    band = [q_start + (q_end - q_start) / (n_q-1) * i for i in range(n_q)]
    return band

# Electronic Modes
kGridDirs = glob("k_grid_conv/*")
kGridDirs.sort(key=lambda s: float(s.split("_")[-1] ) )
kpoints = { "L"      : np.array([0.500, 0.500, 0.500]),
            "\Gamma" : np.array([0.000, 0.000, 0.000]),
            "X"      : np.array([0.000, 0.500, 0.500]),
            "W"      : np.array([0.250, 0.500, 0.750]),
            "K"      : np.array([0.375, 0.375, 0.750])
          }

elecDs = [ DOSFingerprint(True, False, spectraFiles=[s + '/KS_DOS_total.dat'], minE=-15.0, maxE=10.0, nbins=100) for s in kGridDirs ]
elecBs = [ BandStructureFingerprint(True, True, kpoints=kpoints, spectraFiles=[s + '/band1001.out', s + '/band1002.out', s + '/band1003.out', s + '/band1004.out'], minE=-1800.0, maxE=20.0, nbins=20 ) for s in kGridDirs ]

print("DOSFingerprint Convergence")
for ii in range(len(elecDs)):
    print(kGridDirs[ii].split("_")[-1], elecDs[-1].scalar_product(elecDs[ii]))
print("BandStructureFingerprint Convergence")
for ii in range(len(elecBs)):
    print(kGridDirs[ii].split("_")[-1], elecBs[-1].scalar_product(elecBs[ii]))

# Phonons
# Set up the phonopy objects
atoms = read_aims('../phonons/geometry.in')
vol = atoms.get_volume()
smatrix = np.array([[-1,  1,  1],
                    [ 1, -1,  1],
                    [ 1,  1, -1]])

# A series of super cell matrices
smatrices = [ a * smatrix for a in range(1,3)]
phononCalcs = [ ph.preprocess(atoms, smatrix) for smatrix in smatrices ]
workdirs = [Path('./Si_{}{}{}_{}{}{}_{}{}{}_{:.3f}'.format(*smatrix.flatten(), vol)) for smatrix in smatrices ]
for direc in workdirs:
    direc.mkdir(exist_ok=True)

# Calculate the Forces
for ii, cell in enumerate(phononCalcs):
    for jj, cell in enumerate(phononCalcs[ii][2]):
        folder_with_disp = workdirs[ii] / f'disp-{jj:03d}'
        folder_with_disp.mkdir(parents=True, exist_ok=True)
        print(str(folder_with_disp) + ' created')
        write(str(folder_with_disp / 'geometry.in'), phononCalcs[ii][2], 'aims', scaled=True)
        try: (folder_with_disp / 'control.in').unlink()
        except: pass
        shutil.copy('control.in', str(folder_with_disp / 'control.in'))

# Define q_path and q_points
q_points = { '\\Gamma': np.array([0.0, 0.00, 0.00]),
             '\\Delta': np.array([0.25, 0., 0.25]),
                   'X': np.array([0.5, 0.00, 0.50]),
                   'W': np.array([0.5, 0.25, 0.75]),
                   'K': np.array([0.375, 0.375, 0.75]),
            '\\Lambda': np.array([0.25, 0.25, 0.25]),
                   'L': np.array([0.5, 0.5, 0.5])
             }

# This is a possible path through the Brillouin zone
q_path = ["\\Gamma", "\\Delta", "X", "W", "K", "\\Gamma", "\\Lambda", "L"]

fpList = []
# Collect the forces
for ii in range(len(workdirs)):
    force_sets = []
    disps = sorted(workdirs[ii].glob('disp-???'))
    phonon = phononCalcs[ii][0]
    for wd in disps:
        try:
            forces = read_aims_output(str(wd / 'aims.out'))[0].get_forces()
        except FileNotFoundError:
            exit(f'Please calculate the forces in {wd} in order to proceed.')
        force_sets.append(forces)
    print(f'.. {len(force_sets)} force(s) have been read from aims output files.')
    phonon.set_forces(force_sets)
    phonon.produce_force_constants()
    bands = []
    for jj, _ in enumerate(q_path[:-1]):
        q_start, q_end = q_points[q_path[jj]], q_points[q_path[jj+1]]
        band = get_band(q_start, q_end)
        bands.append(band)
    phonon.set_band_structure(bands)
    latexify = lambda sym: "$\\mathrm{\\mathsf{" + str(sym)  + "}}$"
    labels = [latexify(sym) for sym in q_path]
    fp = get_phonon_bs_fingerprint_phononpy(phonon, q_points, vectorize=True)
    fpList.append( fp )
    q_mesh = [45, 45, 45]
    phonon.set_mesh(q_mesh)
    # We generate the DOS by calling .set_total_DOS on the phonopy object
    phonon.set_total_DOS(freq_pitch=.1, tetrahedron_method=True)
print("phonon BS scalar product")
for ff in fpList:
    print( scalar_product(ff, fpList[-1], True, True) )

phonon_dos_fp = get_phonon_dos_fingerprint_phononpy(phononCalcs[0][0])
print("Phonon DOS")
print(phonon_dos_fp["DOS"])