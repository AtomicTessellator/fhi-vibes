from pathlib import Path

import numpy as np
from ase.io import read

from vibes.harmonic_analysis.mode_projection import SimpleModeProjection
from vibes.io import parse_force_constants

parent = Path(__file__).parent

fc_file = parent / "infile.forceconstant"
primitive = read(parent / "geometry.in.primitive", format="aims")
supercell = read(parent / "geometry.in.supercell", format="aims")

ref_omegas = np.loadtxt(parent / "ref_omegas.dat")

fc = parse_force_constants(fc_file, primitive=primitive, supercell=supercell)


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
