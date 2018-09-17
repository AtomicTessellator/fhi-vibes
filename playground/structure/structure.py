"""
Contains routines to manipulate structures, such as
 - symmetry refinement
 - size scaling
 - wrapping
 - possibly supercell generation


 Future idea for the format:
   - sort types of atoms according to atomic number
   - sort according to distance from origin

 Future speed considerations:
   - make symmetry dataset an attribute to reduce no. of calls
"""
from ase.atoms import Atoms
import itertools
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin #faster than numpy for scalars
import datetime
from .symmetry import Spacegroup
from . import structure_io
from playground.konstanten.symmetry import symprec

class Cell(Atoms):
    def __init__(self,
                 ase_atoms=None,
                 phonopy_atoms=None,
                 symprec = symprec,
                 **kwargs):

        if ase_atoms:
            # use the atoms object as provided
            pass
        elif phonopy_atoms:
            # convert to ase.atoms.Atoms first
            ase_atoms = Atoms(symbols   = phonopy_atoms.get_chemical_symbols(),
                              cell      = phonopy_atoms.get_cell(),
                              positions = phonopy_atoms.get_positions())
            ase_atoms.set_pbc(True)
        else:
            # create an Atoms object from the argumets provided
            ase_atoms = Atoms(**kwargs)

        # initialize ase Atoms object
        super(Cell, self).__init__(ase_atoms)

        if symprec:
            self.spacegroup  = Spacegroup(self, symprec)
        else:
            self.spacegroup  = None
        #
        self.symprec     = symprec
        self.sysname     = self.get_sysname()
        self.Natoms      = self.get_number_of_atoms()
        self.tags        = [None]
        #
        # Constraints:
        self.constraints_pos = [None for ii in range(self.Natoms)]
        self.constraints_lv  = [None, None, None]
        self.symmetry_block  = []

    @property
    def n_atoms(self):
        return self.get_number_of_atoms()

    @property
    def symbols(self):
        return self.get_chemical_symbols()

    def get_unique_symbols(self):
        return np.unique(self.symbols, return_counts=True)

    def get_sysname(self):
        """ Get name of the system:
        Either the chemical formula, or the chemical formula enriched by spacegroup information"""

        chemical_formula      = self.get_chemical_formula()
        # if there is no spacegroup
        if self.spacegroup is None:
            return chemical_formula

        sg_number             = self.spacegroup.number
        wyckoff_pos           = self.spacegroup.wyckoffs
        sysname               = f'{chemical_formula}_{sg_number}'
        wyck_uniq, wyck_mult  = np.unique(wyckoff_pos, return_counts=1)
        for mult, wyck in zip(wyck_mult, wyck_uniq):
            sysname += f'_{mult}{wyck}'
        return sysname

    def get_conventional_standardized(self):
        return self.spacegroup.get_conventional_standardized()

    def get_primitive_standardized(self):
        return self.spacegroup.get_primitive_standardized()


    def get_nn_dist(self):
        """ returns nearest neighbor distance in primitive """
        try:
            dists = self.get_all_distances(mic=True)
        except ValueError:
            dists = self.get_conventional_standardized().get_all_distances()
        dists = np.triu(dists)  # only upper triag. because of symmetry
        mindist = np.amin( dists[np.nonzero(dists)] ) # smallest non zero dist.
        return mindist

    def constrain_lv(self, lv = None):
        if lv == 0:
            self.constraints_lv[0] = True
        elif lv == 1:
            self.constraints_lv[1] = True
        elif lv == 2:
            self.constraints_lv[2] = True
    #

    def get_string(self, decorated=True, format = 'aims', scaled=True):
        if format == 'aims':
            return structure_io.get_aims_string(self, decorated = decorated, scaled = True)
        #
        else:
            print(f'Structure output format {format} not implemented. Stop.')
            exit()

    def write(self, filename='geometry.in', format='aims', scaled = True):
        with open(filename, 'w') as f:
            f.write(self.get_string(format=format, scaled = True))
            # Write symmetry block if existent:
            if len(self.symmetry_block) > 0:
                f.write('\n#Symmetry block for constrained relaxation:\n')
                for line in self.symmetry_block:
                    f.write(line)

    def inform(self, *args, **kwargs):
       structure_io.inform(self, *args, **kwargs)

    #
    def get_hash(self, short = False):
        from playground.helpers.hash import hash_atoms
        return hash_atoms(self)[0]

    #
    # Soubroutines to find cubic supercells
    def find_cubic_cell(self, Ntarget=None,
                        deviate_fill = 0.1, deviate_N = 0.1,
                        verbose=False, analytical=False):
        from .cubic_supercell import find_cubic_supercell

        Sfin = find_cubic_supercell(self, Ntarget,
                        deviate_fill, deviate_N,
                        verbose, analytical)
        return Sfin

    def make_supercell(self, P):
        """Generate a supercell by applying a general transformation (*P*) to
        the input configuration (*prim*).

        The transformation is described by a 3x3 integer matrix
        `\mathbf{P}`. Specifically, the new cell metric
        `\mathbf{h}` is given in terms of the metric of the input
        configuraton `\mathbf{h}_p` by `\mathbf{P h}_p =
        \mathbf{h}`.

        Internally this function uses the :func:`~ase.build.cut` function.

        Parameters:

        prim: ASE Atoms object
            Input configuration.
        P: 3x3 integer matrix
            Transformation matrix `\mathbf{P}`.

        """

        from ase.build import cut
        atom = Atoms(symbols   = self.get_chemical_symbols(),
                     cell      = self.get_cell(),
                     positions = self.get_positions())
        atom.set_pbc(True)

        return Cell(cut(atom, P[0], P[1], P[2]))

    def refine(self, primitive=True):
        import spglib as spg
        if not hasattr(self, 'symprec'):
            exit('Structure object does not have symprec attribute, but symmetry' +
                  ' refinement was requested. Abort.')

        lattice, scaled_positions, numbers = spg.standardize_cell(
            self, to_primitive=primitive, no_idealize=0, symprec=self.symprec)
        refined_cell = Atoms(cell=lattice, scaled_positions=scaled_positions,
                             numbers=numbers, pbc=True)

        refined_cell.wrap()
        return Cell(refined_cell, self.symprec)


    def change_volume(self, factor):
        scaled = self.get_scaled_positions()
        newcell = []
        fac = factor**(1./3)
        for latvec in self.cell:
            newcell.append(latvec*fac)
        new_structure = Cell(Atoms(cell=newcell, scaled_positions=scaled,
                        numbers=self.numbers, pbc=True))
        return new_structure

    def make_from_cell(self, newcell):
        scaled = self.get_scaled_positions()
        new_structure = Cell(Atoms(cell=newcell, scaled_positions=scaled,
                        numbers=self.numbers, pbc=True))
        return new_structure
