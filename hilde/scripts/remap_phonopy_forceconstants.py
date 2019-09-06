"""read FORCE_CONSTANTS and remap to given supercell"""

import numpy as np

from hilde.helpers import talk
from hilde.io import parse_force_constants


def remap_force_constants(
    fc_filename="FORCE_CONSTANTS",
    uc_filename="geometry.primitive",
    sc_filename="geometry.supercell",
    fortran=True,
    eps=1e-13,
    tol=1e-5,
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

    fc = parse_force_constants(
        fc_file=fc_filename,
        primitive=uc_filename,
        supercell=sc_filename,
        fortran=fortran,
        two_dim=True,
        eps=eps,
        tol=tol,
    )

    msg = f"remapped force constants from {fc_filename}, shape [{fc.shape}]"
    outfile = f"{fc_filename}_remapped"
    np.savetxt(outfile, fc, header=msg)

    talk(f".. remapped force constants written to {outfile}")
