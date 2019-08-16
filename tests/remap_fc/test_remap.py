from pathlib import Path
import numpy as np
from hilde.phonopy.utils import parse_phonopy_force_constants
from phonopy.file_IO import parse_FORCE_CONSTANTS

parent = Path(__file__).parent

folders = ["."]


def _symmetric(mat):
    violation = np.linalg.norm(mat - mat.T)
    assert violation < 1e-12, violation


def _parse(folder, fortran=True):
    fc = parse_phonopy_force_constants(
        parent / folder / "FORCE_CONSTANTS",
        primitive=parent / folder / "geometry.in.primitive",
        supercell=parent / folder / "geometry.in.supercell",
        fortran=fortran,
    )
    return fc


def _test_folder(folder="."):
    fc_fortran = _parse(folder, fortran=True)
    _symmetric(fc_fortran)

    fc_python = _parse(folder, fortran=False)
    _symmetric(fc_python)

    fc_phonopy = parse_FORCE_CONSTANTS(parent / folder / "FORCE_CONSTANTS_reference")
    fc_phonopy = fc_phonopy.swapaxes(1, 2).reshape(2 * (3 * fc_phonopy.shape[1],))

    norm = np.linalg.norm(fc_python - fc_fortran)
    assert norm < 1e-12, (norm, folder)

    # check phonopy
    norm = np.linalg.norm(fc_python - fc_phonopy)
    assert norm < 1e-12, (norm, folder)


def test_folders(folders=folders):
    for folder in folders:
        _test_folder(folder)


if __name__ == "__main__":
    test_folders()
