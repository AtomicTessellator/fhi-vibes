import sys
sys.path.append('../src/')
from MaterialsFingerprints import MaterialsFingerprint
import numpy as np

# Make initial fingerprints
elecD = MaterialsFingerprint(True, False, spectraFiles=['D_Fingerprint/electron/KS_DOS_total.dat'], dE=0.1, minE=0.0, maxE=6.0)
elecB = MaterialsFingerprint(True, True, kpoints={'\Gamma':np.array([0.0,0.0,0.0]), 'J':np.array([0.5,0.0,0.0]), 'K':np.array([0.5,0.5,0.0])}, spectraFiles=['B_Fingerprint/electron/band1001.out', 'B_Fingerprint/electron/band1002.out'], minE=-14.0, maxE=14.0)

phononD = MaterialsFingerprint(False, False, spectraFiles=['D_Fingerprint/phonon/total_dos.dat'], dE=0.1, minE=0.0, maxE=6.0)
phononB = MaterialsFingerprint(False, True, kpoints={ '\Gamma': np.array([0.0, 0.00, 0.00]), 'X': np.array([0.0, 0.50, 0.0]), 'M': np.array([0.5, 0.5, 0.0]), 'R': np.array([0.5, 0.5, 0.5]) }, spectraYAML='B_Fingerprint/phonon/band.yaml')

print(elecD.scalarProduct(phononD))
