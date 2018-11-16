""" Helpers for working with kpoint densites """

import numpy as np


def d2k(atoms, kptdensity=3.5, even=True):
    """[ase.calculators.calculator.kptdensity2monkhorstpack]
    Convert k-point density to Monkhorst-Pack grid size.

    atoms: Atoms object
        Contains unit cell and information about boundary conditions.
    kptdensity: float
        Required k-point density.  Default value is 3.5 point per Ang^-1.
    even: bool
        Round up to even numbers.
    """

    recipcell = atoms.get_reciprocal_cell()
    kpts = []
    for i in range(3):
        if atoms.pbc[i]:
            k = 2 * np.pi * np.sqrt((recipcell[i] ** 2).sum()) * float(kptdensity)
            if even:
                kpts.append(2 * int(np.ceil(k / 2)))
            else:
                kpts.append(int(np.ceil(k)))
        else:
            kpts.append(1)
    return kpts


def k2d(atoms, k_grid=[2, 2, 2]):
    """Generate the kpoint density in each direction from given k_grid.

    Args:
        atoms (Atoms): Atoms object of interest.
        k_grid (list): k_grid that was used.

    Returns:
        np.ndarray: The density of kpoints in each direction.
                    Use result.mean() to compute average kpoint density.
    """

    recipcell = atoms.get_reciprocal_cell()
    densities = k_grid / (2 * np.pi * np.sqrt((recipcell ** 2).sum(axis=1)))
    return np.array(densities)


def update_k_grid(atoms, calc, kptdensity, even=True):
    """ update the k_grid in calc with the respective density """

    k_grid = d2k(atoms, kptdensity, even)

    if "k_grid" in calc.parameters:
        calc.parameters["k_grid"] = k_grid
