import sys
from hilde.MatRefDatabase.MaterialsFingerprints import MaterialsFingerprint, DOSFingerprint, BandStructureFingerprint
import numpy as np
from glob import glob

kGridDirs = glob("k_grid_conv/*")
kGridDirs.sort(key=lambda s: float(s.split("_")[-1] ) )
kpoints = { "L"      : np.array([0.500, 0.500, 0.500]),
		    "\Gamma" : np.array([0.000, 0.000, 0.000]),
		    "X"      : np.array([0.000, 0.500, 0.500]),
		    "W"      : np.array([0.250, 0.500, 0.750]),
		    "K"      : np.array([0.375, 0.375, 0.750])
		  }
elecDs = [ DOSFingerprint(True, False, spectraFiles=[s + '/KS_DOS_total.dat'], minE=-2.00e3, maxE=2.0e1, nbins=301) for s in kGridDirs ]
elecBs = [ BandStructureFingerprint(True, True, kpoints=kpoints, spectraFiles=[s + '/band1001.out', s + '/band1002.out', s + '/band1003.out', s + '/band1004.out'], minE=-2.0e3, maxE=2.0e1, nbins=51 ) for s in kGridDirs ]
print("DOSFingerprint")
for ii in range(len(elecDs)):
	print(kGridDirs[ii].split("_")[-1], elecDs[-1].scalarProduct(elecDs[ii]))
print("BandStructureFingerprint")
for ii in range(len(elecBs)):
	print(kGridDirs[ii].split("_")[-1], elecBs[-1].scalarProduct(elecBs[ii]))

# Make initial fingerprints
# elecD = DOSFingerprint(True, False, spectraFiles=['D_Fingerprint/electron/KS_DOS_total.dat'], dE=0.1, minE=0.0, maxE=6.0)
# elecB = BandStructureFingerprint(True, True, kpoints={'\Gamma':np.array([0.0,0.0,0.0]), 'J':np.array([0.5,0.0,0.0]), 'K':np.array([0.5,0.5,0.0])}, spectraFiles=['B_Fingerprint/electron/band1001.out', 'B_Fingerprint/electron/band1002.out'], minE=-14.0, maxE=14.0)

# phononD = DOSFingerprint(False, False, spectraFiles=['D_Fingerprint/phonon/total_dos.dat'], dE=0.1, minE=0.0, maxE=6.0)
# phononB = BandStructureFingerprint(False, True, kpoints={ '\Gamma': np.array([0.0, 0.00, 0.00]), 'X': np.array([0.0, 0.50, 0.0]), 'M': np.array([0.5, 0.5, 0.0]), 'R': np.array([0.5, 0.5, 0.5]) }, spectraYAML='B_Fingerprint/phonon/band.yaml')

