""" Module containing wrapper functions to work with Phonopy """

from hilde.structure import pAtoms


def to_pAtoms(phonopy_atoms, smatrix, symprec=None,
              info=None):
    """ Convert one or several PhonopyAtoms to pAtoms

    Args:
        phonopy_atoms (PhonopyAtoms/list): one or several PhonopyAtoms
        smatrix (ndarray): Supercell matrix
        symprec (float): symmetry precision
        info (str): can be used to enumerate phonpopy_displacments

    Returns:
        pAtoms/list: one or several pAtoms

    """

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
        if info is not None:
            atoms.info[info] = ii
        atoms.info['supercell'] = True
        atoms.info['smatrix'] = list(smatrix.T.flatten())
        out_atoms.append(atoms)

    if isinstance(phonopy_atoms, list):
        return out_atoms
    return out_atoms[0]
