from hilde.MatRefDatabase.MaterialsFingerprints import MaterialsFingerprint
from hilde.MatRefDatabase.DatabaseClass import *
from hilde.MatRefDatabase.databaseRegisterFunctions import *
import sqlite3
from ase.io import read
from ase.calculators.aims import Aims
from hilde.structure import Cell
from hilde.structure import read_aims
import numpy as np

# Make initial fingerprints
elecD = MaterialsFingerprint(True, False, spectraFiles=['D_Fingerprint/electron/KS_DOS_total.dat'])
elecB = MaterialsFingerprint(True, True, kpoints={'\Gamma':np.array([0.0,0.0,0.0]), 'J':np.array([0.5,0.0,0.0]), 'K':np.array([0.5,0.5,0.0])}, spectraFiles=['B_Fingerprint/electron/band1001.out', 'B_Fingerprint/electron/band1002.out'], minE=-14.0, maxE=14.0)

phononD = MaterialsFingerprint(False, False, spectraFiles=['D_Fingerprint/phonon/total_dos.dat'])
phononB = MaterialsFingerprint(False, True, kpoints={ '\Gamma': np.array([0.0, 0.00, 0.00]), 'X': np.array([0.0, 0.50, 0.0]), 'M': np.array([0.5, 0.5, 0.0]), 'R': np.array([0.5, 0.5, 0.5]) }, spectraYAML='B_Fingerprint/phonon/band.yaml')

db = MateirailsDatabase("test_fingerprints.db", [(MaterialsFingerprint, adapt_fingerprint), (Atoms, adapt_atom), (Cell, adapt_atom)], [("fingerprint", convert_fingerprint), ("atom", convert_atom), ("cell", convert_atom)])

# db.addAdapter2Register(MaterialsFingerprint, adapt_fingerprint)
# db.addConverter2Register("fingerprint", convert_fingerprint)

if(len(db.getTableList() ) > 0):
	db.dropTable("fingerprints")

db.createTable( "fingerprints", [("id", "integer"),("elecD", "fingerprint"),("elecB", "fingerprint"), ("phononD", "fingerprint"), ("phononB", "fingerprint")])

db.insertRows("fingerprints", 'id', [1, 2], [ [ ("elecD", elecD), ("elecB", elecB) ], [ ("phononD", phononD), ("phononB", phononB) ] ])
db.updateRow('fingerprints', 'id', [1, 2], [ [ ("elecD", elecD), ("elecB", elecB) ], [ ("phononD", phononD), ("phononB", phononB) ] ])

row = db.selectRows('fingerprints', 'id', [1,2])

testFingerPrintRet = db.selectCell("fingerprints", 'id', 1, "elecD")

db.createTable("atoms", [("id", "integer"), ("atomsObj", "atom"), ("element", "string_list")])

print('read_aims:')
atoms = read_aims('../../test/geometry.in')

atoms.calc = Aims(k_grid=[1,1,1],
                  aims_command='orterun -n 4 /home/knoop/FHIaims/bin/aims.ipi.mpi.x',
                  species_dir='/home/knoop/FHIaims/aimsfiles/species_defaults/light/',
                  xc='pw-lda',
                  output_level='MD_light'
                  )
print(type(atoms))
db.insertRows("atoms", 'id', [1], [[("atomsObj", atoms)]])
atomsCopy = db.selectCell("atoms", 'id', 1, "atomsObj")
print(adapt_atom(atoms))
print(adapt_atom(atomsCopy))
db.close()