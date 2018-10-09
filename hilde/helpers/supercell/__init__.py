import numpy as np
from . import supercell as sc
from ase.build.supercells import make_supercell as ase_make_supercell
from hilde.structure import pAtoms

def find_cubic_cell(cell,
                    target_size=1,
                    deviation=0.2,
                    lower_limit=-2,
                    upper_limit=2,
                    verbose=False):

    smatrix = sc.supercell.find_optimal_cell(cell, np.eye(3),
                                             target_size=target_size,
                                             deviation=deviation,
                                             lower_limit=lower_limit, upper_limit=upper_limit,
                                             verbose=verbose
                                             )
    return smatrix

def make_cubic_supercell(atoms,
                         target_size=100,
                         deviation=0.2,
                         lower_limit=-2,
                         upper_limit=2,
                         verbose=False):

    target_size /= len(atoms)

    smatrix = find_cubic_cell(cell=atoms.cell,
                              target_size=target_size,
                              deviation=deviation,
                              lower_limit=lower_limit,
                              upper_limit=upper_limit,
                              verbose=verbose)

    supercell = make_supercell(atoms, smatrix)
    return supercell, smatrix

def make_supercell(*args, **kwargs):
    """ Wrap the make_supercell() function from ase.build """
    supercell = pAtoms(ase_atoms=ase_make_supercell(*args, **kwargs))
    supercell.tags.append('supercell')
    return supercell
