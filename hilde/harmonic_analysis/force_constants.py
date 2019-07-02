""" Tools for dealing with force constants """

import numpy as np
from hilde.helpers.numerics import clean_matrix
from hilde.helpers.lattice_points import get_lattice_points, map_I_to_iL
from hilde.phonopy.utils import remap_force_constants


def reshape_force_constants(
    primitive, supercell, force_constants, scale_mass=False, lattice_points=None
):
    """ reshape from (N_prim x N_super x 3 x 3) into 3x3 blocks labelled by (i,L)

    Parameters
    ----------
    primitive: ase.atoms.Atoms
        The primitive cell structure
    supercell: ase.atoms.Atoms
        The super cell structure
    force_constants: np.ndarray(shape=(N_prim x N_super x 3 x 3))
        The input force constant matrix
    scale_mass: bool
        If True scale phi by the product of the masses
    lattice_points: np.ndarray
        the lattice points to include

    Returns
    -------
    new_force_constants: np.ndarray(shape=(i,L))
        The remapped force constants
    """
    if len(force_constants.shape) > 2:
        if force_constants.shape[0] != force_constants.shape[1]:
            force_constants = remap_force_constants(
                force_constants,
                primitive,
                supercell,
                new_supercell=None,
                reduce_fc=False,
                two_dim=True,
            )
        else:
            force_constants = force_constants.swapaxes(1, 2).reshape(
                2 * (3 * force_constants.shape[0])
            )

    if lattice_points is None:
        lattice_points, _ = get_lattice_points(primitive.cell, supercell.cell)

    indeces = map_I_to_iL(primitive, supercell, lattice_points=lattice_points)

    n_i = len(primitive)
    n_L = len(lattice_points)

    masses = primitive.get_masses()

    new_force_constants = np.zeros([n_i, n_L, n_i, n_L, 3, 3])

    for n1 in range(len(supercell)):
        for n2 in range(len(supercell)):
            phi = force_constants[3 * n1 : 3 * n1 + 3, 3 * n2 : 3 * n2 + 3]

            i1, L1, i2, L2 = (*indeces[n1], *indeces[n2])

            if scale_mass:
                phi /= np.sqrt(masses[i1] * masses[i2])

            new_force_constants[i1, L1, i2, L2] = clean_matrix(phi)

    return new_force_constants
