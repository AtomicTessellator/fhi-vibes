from pathlib import Path

import numpy as np
from ase.io import read

from vibes.force_constants import reduce_force_constants, remap_force_constants
from vibes.helpers.supercell import map2prim
from vibes.phonopy.utils import parse_phonopy_force_constants


parent = Path(__file__).parent


def test_reduce_fc():
    fc_phonopy = parse_phonopy_force_constants(parent / "FORCE_CONSTANTS")

    kw = {"format": "aims"}
    primitive = read(parent / "geometry.in.primitive", **kw)
    supercell = read(parent / "geometry.in.supercell", **kw)

    kw = {"primitive": primitive, "supercell": supercell}
    fc_remapped = remap_force_constants(fc_phonopy, **kw, symmetrize=False)

    sc2pc = map2prim(**kw)
    fc_reduced = reduce_force_constants(fc_remapped, sc2pc)
    print(fc_reduced.shape)

    diff = np.linalg.norm(fc_phonopy - fc_reduced)

    assert np.allclose(diff, 0), f"diff should be 0, is {diff}"


if __name__ == "__main__":
    test_reduce_fc()
