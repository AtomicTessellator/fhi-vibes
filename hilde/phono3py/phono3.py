"""
A leightweight wrapper for Phono3py
"""

from collections import namedtuple
import numpy as np
from phono3py.phonon3 import Phono3py
from hilde import konstanten as const
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


def prepare_phono3py(atoms,
                     fc2_supercell_matrix,
                     fc3_supercell_matrix,
                     q_mesh=[11, 11, 11],
                     fc2=None,
                     fc3=None,
                     disp=0.03,
                     cutoff_pair_distance=10.,
                     symmetrize_fc3q=False,
                     symprec=1e-5,
                     log_level=2):
    """Prepare a Phono3py object.

    Args:
        atoms (pAtoms): (primitive) unit cell
        fc2_supercell_matrix (ndarray): smatrix for FC2 supercell
        fc3_supercell_matrix (ndarray): smatrix for FC3 supercell
        q_mesh (list): q-point mesh for BZ integrations
        fc2 (ndarray): Pre-computed FC2
        fc3 (ndarray): Pre-computed FC3
        disp (float): Finitie displacemt
        cutoff_pair_distance (float): Cutoff radius for triplet interactions
        symmetrize_fc3q (bool): symmetrize fc in q space
        symprec (float): Symmetry tolerance
        log_level (int): Phono3py log level

    Returns:
        Phono3py(): a phono3py object

    """
    ph_atoms = atoms.to_phonopy_atoms()

    phonon3 = Phono3py(ph_atoms,
                       supercell_matrix=fc3_supercell_matrix,
                       phonon_supercell_matrix=fc2_supercell_matrix,
                       mesh=q_mesh,
                       symprec=symprec,
                       is_symmetry=True,
                       symmetrize_fc3q=symmetrize_fc3q,
                       frequency_factor_to_THz=const.eV_to_THz,
                       log_level=log_level)

    phonon3.generate_displacements(distance=disp,
                                   cutoff_pair_distance=cutoff_pair_distance)

    if fc2 is not None:
        phonon3.set_fc2(fc2)
    if fc3 is not None:
        phonon3.set_fc3(fc3)

    return phonon3


def preprocess(atoms,
               fc2_supercell_matrix,
               fc3_supercell_matrix,
               q_mesh=[11, 11, 11],
               disp=0.03,
               cutoff_pair_distance=10.,
               symprec=1e-5):
    """
    Set up a Phono3py object and generate all the necessary supercells

    Args:
        atoms: atoms object that represents the (primitive) unit cell
        q_mesh: q-point grid for BZ integrations
        supercell_matrix: supercell matrix
        disp: displacement for the finite displacemt

    Returns:
        namedtuple with the phonon3 object, the supercell
        and the supercells_with_displacements as hilde.pAtoms
    """

    phonon3 = prepare_phono3py(atoms,
                               fc2_supercell_matrix=fc2_supercell_matrix,
                               fc3_supercell_matrix=fc3_supercell_matrix,
                               q_mesh=q_mesh,
                               disp=disp,
                               cutoff_pair_distance=cutoff_pair_distance,
                               symprec=symprec)

    # phonpoy supercells
    fc2_supercell = to_pAtoms(phonon3.get_phonon_supercell(),
                              fc2_supercell_matrix,
                              symprec=symprec)

    fc3_supercell = to_pAtoms(phonon3.get_supercell(),
                              fc3_supercell_matrix,
                              symprec=symprec)

    scells = phonon3.get_phonon_supercells_with_displacements()
    fc2_supercells_with_disps = to_pAtoms(scells, fc2_supercell_matrix)

    scells = phonon3.get_supercells_with_displacements()
    fc3_supercells_with_disps = to_pAtoms(scells, fc3_supercell_matrix)

    pp = namedtuple('phono3py_preprocess',
                    ('phonon3 fc2_supercell fc3_supercell' +
                     ' fc2_supercells_with_displacements' +
                     ' fc3_supercells_with_displacements'))

    return pp(phonon3, fc2_supercell, fc3_supercell,
              fc2_supercells_with_disps, fc3_supercells_with_disps)


def get_forces(supercells, supercells_computed):
    """ Return force_sets taking care of supercells that were not computed
    because of cutoff. """

    zero_force = np.zeros([len(supercells[0]), 3])
    force_sets = []
    counter = 0
    for scell in supercells:
        if scell is None:
            force_sets.append(zero_force)
        else:
            force_sets.append(supercells_computed[counter].get_forces())
            counter += 1

    if len(force_sets) != len(supercells):
        print('len(force_sets), len(supercells):',
              len(force_sets), len(supercells))
        raise RuntimeError("Number of computed supercells incorrect.")

    return force_sets

# def get_force_constants(phonon, force_sets=None):
#     """
#     fkdev: is this necessary?
#     Take a Phonopy object and produce force constants from the given forces
#     """
#     n_atoms = phonon.get_supercell().get_number_of_atoms()
#
#     phonon.produce_force_constants(force_sets)
#
#     force_constants = phonon.get_force_constants()
#
#     if force_constants is not None:
#         # convert forces from (N, N, 3, 3) to (3*N, 3*N)
#         force_constants = phonon.get_force_constants().swapaxes(1, 2).reshape(2*(3*n_atoms, ))
#         return force_constants
#     # else
#     warn('Force constants not yet created, please specify force_sets.')
