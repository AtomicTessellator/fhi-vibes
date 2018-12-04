from structure import read_output
from ase.io import read
from helpers.config import ConfigDict
from helpers.hash import hash_atoms_and_calc, atoms2dict
from ase.calculators.aims import Aims

atoms = read('./test/geometry.in', 0, 'aims')
print('ase read:')
print(hash_atoms_and_calc(atoms))

print('read_aims:')
atoms = read_output('./test/geometry.in')
print(hash_atoms_and_calc(atoms))

atoms.calc = Aims(k_grid=[1,1,1],
                  aims_command='orterun -n 4 /home/knoop/FHIaims/bin/aims.ipi.mpi.x',
                  species_dir='/home/knoop/FHIaims/aimsfiles/species_defaults/light/',
                  xc='pw-lda',
                  output_level='MD_light'
                  )

print('after calc attached:')
print(hash_atoms_and_calc(atoms))

atomsdict = atoms2dict(atoms)

print('after ignore list:')
print(hash_atoms_and_calc(atoms,
                 ignore_file='./test/hash_ignore.ini'))

atoms.calc = Aims(k_grid=[1,1,1],
                  aims_command='orterun -n 4 /home/knoop/FHIaims/bin/aims.ipi.mpi.x',
                  species_dir='/home/knoop/FHIaims/aimsfiles/species_defaults/light/',
                  xc='pw-lda'
                  )

print('after calc w/o output_level attached:')
print(hash_atoms_and_calc(atoms))
