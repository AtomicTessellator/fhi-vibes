"""read FORCE_CONSTANTS and remap to given supercell"""

import numpy as np

from hilde.helpers import talk
from hilde.phonopy.utils import parse_phonopy_force_constants


def remap_phonopy_force_constants(
    uc_filename="geometry.primitive",
    sc_filename="geometry.supercell",
    fc_filename="FORCE_CONSTANTS",
    eps=1e-13,
    tol=1e-5,
    format="aims",
):
    """take phonopy FORCE_CONSTANTS and remap to given supercell

    Parameters
    ----------
    uc_filename: str or Path
        The input file for the primitive/unit cell
    sc_filename: str or Path
        The input file for the supercell
    fc_filename: str or Path
        The phonopy forceconstant file to parse
    eps: float
        finite zero
    tol: float
        tolerance to discern pairs
    format: str
        File format for the input geometries
    """

    fc = parse_phonopy_force_constants(
        uc_filename=uc_filename,
        sc_filename=sc_filename,
        fc_filename=fc_filename,
        two_dim=True,
        eps=eps,
        tol=tol,
        format=format,
    )

    msg = f"remapped force constants from {fc_filename}, shape [{fc.shape}]"
    outfile = f"{fc_filename}_remapped"
    np.savetxt(outfile, fc, header=msg)

    talk(f".. remapped force constants written to {outfile}")
