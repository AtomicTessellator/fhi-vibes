""" helpers to deal with structures (as represented by ase.atoms.Atoms) """

import numpy as np
from ase.geometry import cell_to_cellpar, cellpar_to_cell

from vibes.helpers.numerics import clean_matrix


def clean_atoms(input_atoms, align=False, decimals=10):
    """Clean all arrays in the atoms object up to defined tolerance

    uses `ndarray.round(decimals=decimals) on the arrays"""

    atoms = input_atoms.copy()

    atoms.cell = np.asarray(atoms.cell).round(decimals=decimals)
    spos = atoms.get_scaled_positions(wrap=False)
    atoms.set_scaled_positions(np.asarray(spos).round(decimals=decimals))

    return atoms
