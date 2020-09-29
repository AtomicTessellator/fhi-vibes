"""obtain standardised stresses from atoms"""

import numpy as np
from ase.calculators.calculator import PropertyNotImplementedError
from ase.constraints import full_3x3_to_voigt_6_stress, voigt_6_to_full_3x3_stress

from vibes.helpers import talk, warn


def get_stresses(atoms):
    """Obtain intensive Nx3x3 stresses"""

    try:
        stress = atoms.get_stress(voigt=False)
    except PropertyNotImplementedError:
        raise RuntimeError("Cannot obtain `stresses` if `stress` is not available.")

    stresses = enforce_3x3(atoms.get_stresses())
    stresses = enforce_intensive(stress, stresses, atoms.get_volume())

    return stresses


def enforce_3x3(stresses):
    """Make stresses into Nx3x3 regardless of input form"""

    if is_voigt(stresses):
        return voigt_6_to_full_3x3_stress(stresses)
    else:
        return stresses


def enforce_intensive(stress, stresses, volume):
    """Ensure that sum(stresses) = stress"""

    summed_stresses = stresses.sum(axis=0)

    if np.allclose(stress, summed_stresses):
        return stresses

    else:
        summed_stresses /= volume
        diff = np.linalg.norm(stress - summed_stresses)

        if diff > 1e-7:
            msg = f"Stress and stresses differ by {diff}"
            warn(msg, level=1)
            talk(full_3x3_to_voigt_6_stress(stress), prefix="stress")
            talk(full_3x3_to_voigt_6_stress(summed_stresses), prefix="stresses")

        return stresses / volume


def is_voigt(stresses):
    """True if stresses is n_atoms x 6

    Also asserts that the array has the
    expected shape.
    """

    if stresses.shape[1] == 6:
        assert len(stresses.shape) == 2

        return True
    else:
        assert stresses.shape[1::] == (3, 3)

        return False
