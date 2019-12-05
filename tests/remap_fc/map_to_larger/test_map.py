from pathlib import Path

import numpy as np
from ase.io import read as _read
from vibes.io import parse_force_constants
from vibes.phonopy.utils import remap_force_constants

parent = Path(__file__).parent

pairs = (
    ("geometry.in.primitive", "FORCE_CONSTANTS_remapped"),
    ("geometry.in.supercell_4", "FORCE_CONSTANTS_remapped_4"),
    ("geometry.in.supercell_8", "FORCE_CONSTANTS_remapped_8"),
)


def read(*args, **kwargs):
    kwargs.update({"format": "aims"})
    return _read(*args, **kwargs)


def _test_pair(pair):
    primitive = parent / "geometry.in.primitive"
    supercell = parent / "geometry.in.supercell"
    new_supercell = parent / pair[0]

    fc_file = parent / "FORCE_CONSTANTS"
    ref_fc_file = parent / pair[1]

    uc = read(primitive, format=format)
    sc = read(supercell, format=format)
    nsc = read(new_supercell, format=format)

    kwargs = {"primitive": uc, "supercell": sc, "fortran": True}

    fc = parse_force_constants(fc_file, two_dim=False, **kwargs)

    kwargs.update({"new_supercell": nsc, "two_dim": True})

    fc = remap_force_constants(fc, **kwargs)

    fc_ref = np.loadtxt(ref_fc_file)

    diff = np.linalg.norm(fc - fc_ref)
    assert diff < 1e-15, diff


def test(pairs=pairs):

    for pair in pairs:
        _test_pair(pair)


if __name__ == "__main__":
    test()
