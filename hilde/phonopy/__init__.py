""" Module containing wrapper functions to work with Phonopy """

from hilde.structure import pAtoms


def to_pAtoms(phonopy_atoms, smatrix, symprec=None):
    """ Convert one or several PhonopyAtoms to pAtoms

    Args:
        phonopy_atoms (PhonopyAtoms/list): one or several PhonopyAtoms
        smatrix (ndarray): Supercell matrix
        symprec (float): symmetry precision

    Returns:
        pAtoms/list: one or several pAtoms

    """

    if isinstance(phonopy_atoms, list):
        latoms = phonopy_atoms
    else:
        latoms = [phonopy_atoms]

    out_atoms = []
    for atoms in latoms:
        # Check if the atoms does exist (important when working with cutoffs)
        if atoms is None:
            out_atoms.append(None)
            continue

        tags = ['supercell', ('smatrix', list(smatrix.T.flatten()))]
        out_atoms.append(
            pAtoms(phonopy_atoms=atoms,
                   symprec=symprec,
                   tags=tags)
            )

    if isinstance(phonopy_atoms, list):
        return out_atoms
    return out_atoms
