import sys;
sys.path.append('../src/');
from MaterialsFingerprints import *
from DatabaseClass import *
import sqlite3

# Make initial fingerprints
elecD = MaterialsFingerprint(True, False, spectraFiles=['D_Fingerprint/electron/KS_DOS_total.dat'])
elecB = MaterialsFingerprint(True, True, kpoints={'\Gamma':np.array([0.0,0.0,0.0]), 'J':np.array([0.5,0.0,0.0]), 'K':np.array([0.5,0.5,0.0])}, spectraFiles=['B_Fingerprint/electron/band1001.out', 'B_Fingerprint/electron/band1002.out'], minE=-14.0, maxE=14.0)

phononD = MaterialsFingerprint(False, False, spectraFiles=['D_Fingerprint/phonon/total_dos.dat'])
phononB = MaterialsFingerprint(False, True, kpoints={ '\Gamma': np.array([0.0, 0.00, 0.00]), '\Delta': np.array([0.25, 0., 0.25]), 'X': np.array([0.5, 0.00, 0.50]), 'W': np.array([0.5, 0.25, 0.75]), 'K': np.array([0.375, 0.375, 0.75]), '\Lambda': np.array([0.25, 0.25, 0.25]), 'L': np.array([0.5, 0.5, 0.5]) }, spectraYAML='B_Fingerprint/phonon/band.yaml')

db = MateirailsDatabase("test_fingerprints.db", [(MaterialsFingerprint, adapt_fingerprint)], [("fingerprint", convert_fingerprint)])

# db.addAdapter2Register(MaterialsFingerprint, adapt_fingerprint)
# db.addConverter2Register("fingerprint", convert_fingerprint)

if(len(db.getTableList() ) > 0):
	db.dropTable("fingerprints")

db.createTable( "fingerprints", [("id", "integer"),("elecD", "fingerprint"),("elecB", "fingerprint"), ("phononD", "fingerprint"), ("phononB", "fingerprint")])

db.insertRows("fingerprints", 'id', [1, 2], [ [ ("elecD", elecD), ("elecB", elecB) ], [ ("phononD", phononD), ("phononB", phononB) ] ])
db.updateRow('fingerprints', 'id', [1, 2], [ [ ("elecD", elecD), ("elecB", elecB) ], [ ("phononD", phononD), ("phononB", phononB) ] ])

row = db.selectRows('fingerprints', 'id', [1,2])

testFingerPrintRet = db.selectCell("fingerprints", 'id', 1, "elecD")
print(type(testFingerPrintRet))
print(testFingerPrintRet.fingerprint["DOS"] - elecD.fingerprint["DOS"])

# print(type(testFingerPrintRet))
# print(testFingerPrintRet.fingerprint['DOS'])
db.close()