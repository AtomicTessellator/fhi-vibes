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
import datetime
from math import sqrt, pi, cos, sin #faster than numpy for scalars
import itertools
from copy import copy, deepcopy
import numpy as np
from ase.atoms import Atoms
from ase.db.row import atoms2dict, AtomsRow
from .symmetry import Spacegroup
from .misc import get_sysname
from . import io
from hilde.konstanten.symmetry import symprec
from hilde.helpers.maths import clean_matrix
from hilde.konstanten.numerics import loose_tol

def patoms2dict(atoms):
    '''
    Converts a pAtoms object into a dict
    Args:
        atoms: pAtoms or Atoms object
            The pAtoms or Atoms object to be converted into a dictionary
    Returns: atoms_dict (dict)
        The dictionary of atoms
    '''
    atoms_dict = atoms2dict(atoms)
    atoms_dict['info'] = atoms.info
    return atoms_dict

def dict2patoms(atoms_dict):
    '''
    Converts a dict into a pAtoms object
    Args:
        atoms_dict: dict
            A dictionary representing the pAtoms object
    Returns: pAtoms
        The corresponding pAtoms object
    '''
    atoms = pAtoms(AtomsRow(atoms_dict).toatoms(attach_calculator=True))
    if "info" in atoms_dict:
        atoms.info = atoms_dict['info']
    return atoms

class pAtoms(Atoms):
    def __init__(self,
                 ase_atoms=None,
                 phonopy_atoms=None,
                 symprec = symprec,
                 tags = [],
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
        super().__init__(ase_atoms)

        # clean lattice
        self.cell = clean_matrix(self.cell, eps=loose_tol)

        if symprec and len(self) <= 1000:
            self.spacegroup  = Spacegroup(self, symprec)
        elif symprec and len(self) > 1000:
            self.spacegroup = None
            print('**Warning: spacegroup for pAtoms with more than 1000 NOT ' +
                  'computed on default.')
        else:
            self.spacegroup  = None
        #
        self.symprec     = symprec
        self.Natoms      = self.get_number_of_atoms()
        self.tags        = tags
        #
        # Constraints:
        self.constraints_pos = [None for ii in range(self.Natoms)]
        self.constraints_lv  = [None, None, None]
        self.symmetry_block  = []
        self.calc_id = None
    def copy(self):
        new_atoms = super().copy()
        new_atoms.spacegroup = copy(self.spacegroup)
        return new_atoms

    @property
    def n_atoms(self):
        return self.get_number_of_atoms()

    @property
    def sysname(self):
        return get_sysname(self)

    def to_phonopy_atoms(self, wrap=False):
        from .convert import to_phonopy_atoms
        return to_phonopy_atoms(self,wrap=False)

    def to_spglib_cell(self):
        from .convert import to_spglib_cell
        return to_spglib_cell(self)

    def get_unique_symbols(self):
        return np.unique(self.symbols, return_counts=True)

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

    def get_string(self, decorated=True, format = 'aims', scaled=None, wrap=True):
        if format == 'aims':
            return io.get_aims_string(self,
                                      decorated=decorated,
                                      scaled=scaled,
                                      wrap=wrap)
        #
        else:
            print(f'Structure output format {format} not implemented. Stop.')
            exit()

    def write(self, filename='geometry.in', format='aims', scaled=None, wrap=True):
        with open(filename, 'w') as f:
            f.write(self.get_string(format=format, scaled=scaled, wrap=wrap))
            # Write symmetry block if existent:
            if len(self.symmetry_block) > 0:
                f.write('\n#Symmetry block for constrained relaxation:\n')
                for line in self.symmetry_block:
                    f.write(line)

    def inform(self, *args, **kwargs):
       io.inform(self, *args, **kwargs)

    #
    def get_hash(self, short = False):
        from hilde.helpers.hash import hash_atoms
        return hash_atoms(self)[0]

    #
    # Soubroutines to find cubic supercells
    def find_cubic_cell(self, Ntarget=100, deviation=0.2,
                        lower_limit=-2, upper_limit=2,
                        verbose=False):
        from hilde.helpers.supercell import find_cubic_cell as find_cc

        target_size = Ntarget / self.n_atoms

        smatrix = find_cc(self.cell,
                          target_size=target_size,
                          deviation=deviation,
                          lower_limit=lower_limit, upper_limit=upper_limit,
                          verbose=verbose)
        return smatrix

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

        return pAtoms(cut(atom, P[0], P[1], P[2]))

    def refine(self, primitive=True):
        import spglib as spg
        if not hasattr(self, 'symprec'):
            exit('Structure object does not have symprec attribute, ' +
                 'but symmetry refinement was requested. Abort.')

        lattice, scaled_positions, numbers = spg.standardize_cell(
            self, to_primitive=primitive, no_idealize=0, symprec=self.symprec)
        refined_cell = Atoms(cell=lattice, scaled_positions=scaled_positions,
                             numbers=numbers, pbc=True)

        refined_cell.wrap()
        return pAtoms(refined_cell, self.symprec)


    def change_volume(self, factor):
        scaled = self.get_scaled_positions()
        newcell = []
        fac = factor**(1./3)
        for latvec in self.cell:
            newcell.append(latvec*fac)
        new_structure = pAtoms(Atoms(cell=newcell, scaled_positions=scaled,
                        numbers=self.numbers, pbc=True))
        return new_structure

    def make_from_cell(self, newcell):
        scaled = self.get_scaled_positions()
        new_structure = pAtoms(Atoms(cell=newcell, scaled_positions=scaled,
                        numbers=self.numbers, pbc=True))
        return new_structure

    def set_calc_id(self, id_val):
        self.calc_id = id_val
