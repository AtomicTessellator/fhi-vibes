""" helpers to deal with structures (as represented by ase.atoms.Atoms) """

import numpy as np
from ase.geometry import cell_to_cellpar, cellpar_to_cell
from hilde.helpers.numerics import clean_matrix


def clean_atoms(input_atoms, align=False, decimals=10):
    """Clean all arrays in the atoms object up to defined tolerance

    uses `ndarray.round(decimals=decimals) on the arrays"""

    atoms = input_atoms.copy()

    atoms.cell = np.asarray(atoms.cell).round(decimals=decimals)
    atoms.positions = np.asarray(atoms.positions).round(decimals=decimals)

    return atoms
