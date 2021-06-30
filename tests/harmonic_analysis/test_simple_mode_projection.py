from pathlib import Path

import numpy as np
from ase.io import read

from vibes.io import parse_force_constants
from vibes.simple_mode_projection import SimpleModeProjection


parent = Path(__file__).parent

fc_file = parent / "FORCE_CONSTANTS_tdep"  # "infile.forceconstant"
primitive = read(parent / "geometry.in.primitive", format="aims")
supercell = read(parent / "geometry.in.supercell", format="aims")

ref_omegas = np.loadtxt(parent / "ref_omegas.dat")

fc = parse_force_constants(
    fc_file, primitive=primitive, supercell=supercell, two_dim=True
)


def test_instantiation():
    proj = SimpleModeProjection(supercell, fc)

    assert proj.atoms == supercell
    assert np.allclose(proj.force_constants, fc)


def test_omegas():
    proj = SimpleModeProjection(supercell, fc)

    omegas = proj.omegas

    assert np.allclose(omegas[3:], ref_omegas[3:])


if __name__ == "__main__":
    test_instantiation()
    test_omegas()
