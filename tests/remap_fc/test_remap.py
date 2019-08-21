from pathlib import Path
import numpy as np
from hilde.phonopy.utils import parse_phonopy_force_constants, remap_force_constants
from phonopy.file_IO import parse_FORCE_CONSTANTS

from ase.io import read

parent = Path(__file__).parent

folders = ["."]  # , "CdYF3"]


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


def _parse_phonopy(file):
    fc_out = parse_FORCE_CONSTANTS(file)
    fc_out = fc_out.swapaxes(1, 2).reshape(2 * (3 * fc_out.shape[1],))

    # symmetrize
    violation = np.linalg.norm(fc_out - fc_out.T)
    if violation > 1e-5:
        msg = f"**Phonopy force constants are not symmetric by {violation:.2e}."
        msg += " Symmetrize."
        print(msg, flush=True)
        fc_out = 0.5 * (fc_out + fc_out.T)

    _symmetric(fc_out)

    return fc_out

def _remap_fc(folder, fortran):
    force_constants = parse_FORCE_CONSTANTS(parent / folder / "FORCE_CONSTANTS")
    primitive = read(parent / folder / "geometry.in.primitive", format="aims")
    supercell = read(parent / folder / "geometry.in.supercell", format="aims")
    new_supercell = read(parent / folder / "geometry.in.new_supercell", format="aims")

    return remap_force_constants(
        force_constants,
        primitive,
        supercell,
        new_supercell,
        reduce_fc=True,
        fortran=fortran,
    )

def _test_folder(folder="."):
    fc_fortran = _parse(folder, fortran=True)
    _symmetric(fc_fortran)

    fc_python = _parse(folder, fortran=False)
    _symmetric(fc_python)

    fc_phonopy = _parse_phonopy(parent / folder / "FORCE_CONSTANTS_reference")

    norm = np.linalg.norm(fc_python - fc_fortran)
    assert norm < 1e-12, (norm, folder)

    # check phonopy
    norm = np.linalg.norm(fc_python - fc_phonopy)
    assert norm < 1e-12, (norm, folder)

    fc_python_remap = _remap_fc(folder, False)
    fc_fortran_remap = _remap_fc(folder, True)

    norm = np.linalg.norm(fc_python_remap - fc_fortran_remap)
    assert norm < 1e-12


def test_folders(folders=folders):
    for folder in folders:
        _test_folder(folder)


if __name__ == "__main__":
    test_folders()
