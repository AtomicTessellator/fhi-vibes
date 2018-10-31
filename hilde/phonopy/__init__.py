""" Module containing wrapper functions to work with Phonopy """

import numpy as np
from phonopy.structure.atoms import PhonopyAtoms
from hilde.structure import pAtoms


displacement_id_str = 'displacement_id'


def to_phonopy_atoms(atoms):
    phonopy_atoms= PhonopyAtoms(
        symbols   = atoms.get_chemical_symbols(),
        cell      = atoms.get_cell(),
        masses    = atoms.get_masses(),
        positions = atoms.get_positions(wrap=True))
    return phonopy_atoms


def to_pAtoms(phonopy_atoms, smatrix, symprec=None):
    """ Convert one or several PhonopyAtoms to pAtoms

    Args:
        phonopy_atoms (PhonopyAtoms/list): one or several PhonopyAtoms
        smatrix (ndarray): Supercell matrix
        symprec (float): symmetry precision

    Returns:
        pAtoms/list: one or several pAtoms

    """

    smatrix = np.asarray(smatrix)

    if isinstance(phonopy_atoms, list):
        latoms = phonopy_atoms
    else:
        latoms = [phonopy_atoms]

    out_atoms = []
    for ii, atoms in enumerate(latoms):
        # Check if the atoms does exist (important when working with cutoffs)
        if atoms is None:
            out_atoms.append(None)
            continue

        tags = ['supercell', ('smatrix', list(smatrix.T.flatten()))]
        atoms = pAtoms(phonopy_atoms=atoms, symprec=symprec, tags=tags)
        atoms.info['supercell'] = True
        atoms.info['smatrix'] = list(smatrix.T.flatten())
        out_atoms.append(atoms)

    if isinstance(phonopy_atoms, list):
        return out_atoms
    return out_atoms[0]

def enumerate_displacements(cells, info_str=displacement_id_str):
    """ Assign a displacemt id to every atoms obect in cells.

    Args:
        cells (list): atoms objects created by, e.g., phonopy
        info_str (str): how to name the child

    Returns:
        list: cells with id attached to atoms.info (inplace)

    """
    for nn, scell in enumerate(cells):
        if scell is None:
            continue
        scell.info[info_str] = nn
