import os
from ase.calculators.aims import Aims as ASEAims

class Aims(ASEAims):
	def __init__(self, *args, **kwargs):
		super(Aims, self).__init__(*args, **kwargs)

	def write_input(self, atoms, properties=None, system_changes=None,
                    ghosts=None, scaled=False):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        have_lattice_vectors = atoms.pbc.any()
        have_k_grid = ('k_grid' in self.parameters or
                       'kpts' in self.parameters)
        if have_lattice_vectors and not have_k_grid:
            raise RuntimeError('Found lattice vectors but no k-grid!')
        if not have_lattice_vectors and have_k_grid:
            raise RuntimeError('Found k-grid but no lattice vectors!')
        # write_aims(os.path.join(self.directory, 'geometry.in'), atoms, scaled, ghosts)
        atoms.write(os.path.join(self.directory, 'geometry.in'))
        self.write_control(atoms, os.path.join(self.directory, 'control.in'))
        self.write_species(atoms, os.path.join(self.directory, 'control.in'))
        self.parameters.write(os.path.join(self.directory, 'parameters.ase'))
