"""a wrapper for TDEP"""
import numpy as np
from ase.io import read

from vibes.force_constants import ForceConstants
from vibes.helpers import Timer


def parse_tdep_remapped_forceconstant(
    file="outfile.forceconstant_remapped",
    primitive="infile.ucposcar",
    supercell="infile.ssposcar",
    return_array=False,
    format="vasp",
):
    """parse the remapped forceconstants from TDEP and return as ForceConstants object

    Args:
        fc_file: remaped force constants
        primitive: primitive cell file
        supercell: supercell fil

    """
    timer = Timer()

    uc = read(primitive, format=format)
    sc = read(supercell, format=format)

    print(f"Parse force constants from\n  {file}")

    with open(file) as fo:
        n_atoms = int(next(fo).split()[0])
        cutoff = float(next(fo).split()[0])

        assert n_atoms == len(sc), "check supercell size"

        print(f".. Number of atoms:   {n_atoms}")
        print(rf".. Real space cutoff: {cutoff:.3f} \AA")

        # Not yet clear how many lattice points / force constants we will get
        fc = np.zeros([n_atoms, n_atoms, 3, 3])

        for i1 in range(n_atoms):
            n_neighbors = int(next(fo).split()[0])
            for _ in range(n_neighbors):
                # neighbour index
                i2 = int(next(fo).split()[0]) - 1

                # lattice vector
                _ = np.array(next(fo).split(), dtype=float)

                # the force constant matrix for pair (i1, i2)
                phi = np.array([next(fo).split() for _ in range(3)], dtype=float)

                # print(i1, i2, phi)

                fc[i1, i2, :, :] = phi

    fcs = ForceConstants(fc, primitive=uc, supercell=sc)
    timer()

    if return_array:
        return fcs.fc_phonopy

    return fcs
