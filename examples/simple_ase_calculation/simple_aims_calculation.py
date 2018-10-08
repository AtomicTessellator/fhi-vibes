from ase.io import read
from ase.calculators.aims import Aims
from hilde.helpers.paths import cwd
from hilde.tasks.calculate import calculate

atoms = read('si.in', 0, 'aims')

command = 'orterun -n 4 aims.x'
species_dir = '/home/knoop/FHIaims/aimsfiles/species_defaults/light'

aims_settings = {
    'command': command,
    'species_dir': species_dir,
    'output_level': 'MD_light',
    'relativistic': 'atomic_zora scalar',
    'xc': 'pw-lda',
    'k_grid': 3 * [2]
}

calc = Aims(**aims_settings)


# option 1:
atoms.calc = calc
with cwd('tmp', mkdir=True):
    atoms.calc.calculate()
print(atoms.get_total_energy())

# option 2:
err = calculate(atoms, calc, 'tmp_hilde')
if not err:
    print(atoms.get_total_energy())
