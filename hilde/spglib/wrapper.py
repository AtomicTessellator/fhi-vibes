""" a light wrapper for spglib """

import numpy as np
import spglib as spg
from hilde.structure.convert import to_spglib_cell
from hilde.konstanten.symmetry import symprec as default_symprec
from hilde.helpers.config import AttributeDict


def get_symmetry_dataset(atoms, symprec=default_symprec):
    """ return the spglib symmetry dataset """

    dataset = spg.get_symmetry_dataset(to_spglib_cell(atoms), symprec=symprec)

    return AttributeDict(dataset)


def map_unique_to_atoms(atoms, symprec=default_symprec):
    """ map each symmetry unique atom to other atoms as used by phonopy PDOS """

    ds = get_symmetry_dataset(atoms, symprec=symprec)

    uniques = np.unique(ds.equivalent_atoms)

    mapping = [[] for _ in range(len(uniques))]

    for ii, index in enumerate(ds.equivalent_atoms):
        for jj, unique in enumerate(uniques):
            if index == unique:
                mapping[jj].append(ii)

    return mapping
