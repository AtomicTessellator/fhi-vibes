""" helpers to deal with structures (as represented by ase.Atoms) """

import numpy as np
from ase.geometry import cell_to_cellpar, cellpar_to_cell
from hilde.helpers.maths import clean_matrix
from hilde.konstanten.numerics import eps

# no ase
def clean_atoms(input_atoms):
    """Objective: Put position of atom 0 to origin, 
    align 1. lattice vector with x axis, 2. to xy plane
    rotation: change lattice via rotations, else: change via cellpars
    """

    atoms = input_atoms.copy()

    atoms.positions -= atoms.positions[0]

    atoms.wrap()

    scaled_pos = atoms.get_scaled_positions()

    old_lattice = atoms.cell

    cell_params = cell_to_cellpar(old_lattice)

    new_lattice = clean_matrix(cellpar_to_cell(cell_params))

    # Sanity check
    vol0 = np.linalg.det(old_lattice)
    vol1 = np.linalg.det(new_lattice)
    assert abs(vol0 - vol1) < eps, (vol0, vol1)

    atoms.cell = new_lattice
    atoms.set_scaled_positions(scaled_pos)

    atoms.positions = clean_matrix(atoms.positions)

    return atoms
