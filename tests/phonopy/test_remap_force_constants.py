"""test for vibes.phonopy.utils.parse_phonopy_force_constants"""

from pathlib import Path

import numpy as np
from ase.io import read

from vibes.harmonic_analysis.dynamical_matrix import get_frequencies
from vibes.phonopy.utils import parse_phonopy_force_constants

parent = Path(__file__).parent
assets = parent / "assets_remap"
uc_file = assets / "geometry.in.primitive"
sc_file = assets / "geometry.in.supercell"
fc_file = assets / "FORCE_CONSTANTS"

frequencies = np.loadtxt(assets / "frequencies.dat")


def test_remap():
    """Test parsing and remapping force constants"""
    atoms = read(sc_file, format="aims")

    fc = parse_phonopy_force_constants(
        fc_file=fc_file,
        primitive=uc_file,
        supercell=sc_file,
        two_dim=True,
        format="aims",
    )

    freqs = get_frequencies(fc, masses=atoms.get_masses())

    assert np.linalg.norm(freqs - frequencies) < 1e-10, np.linalg.norm(
        freqs - frequencies
    )


if __name__ == "__main__":
    test_remap()
