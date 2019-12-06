"""test for vibes.phonopy.utils.parse_phonopy_force_constants"""
from pathlib import Path
import numpy as np

from ase.io import read
from vibes.phonopy.utils import parse_phonopy_force_constants
from vibes.harmonic_analysis.dynamical_matrix import get_frequencies

parent = Path(__file__).parent
assets = parent / "assets_remap"
uc_filename = assets / "geometry.in.primitive"
sc_filename = assets / "geometry.in.supercell"
fc_filename = assets / "FORCE_CONSTANTS"

frequencies = np.loadtxt(assets / "frequencies.dat")


def test_remap():
    """test parsing and remapping force constants"""

    atoms = read(sc_filename, format="aims")

    fc = parse_phonopy_force_constants(
        fc_filename=fc_filename,
        primitive=uc_filename,
        supercell=sc_filename,
        two_dim=True,
        format="aims",
    )

    freqs = get_frequencies(fc, masses=atoms.get_masses())

    assert np.linalg.norm(freqs - frequencies) < 1e-10, np.linalg.norm(freqs - frequencies)


if __name__ == "__main__":
    test_remap()
