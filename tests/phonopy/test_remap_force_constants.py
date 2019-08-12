from pathlib import Path
import numpy as np
from ase.io import read
from hilde.phonopy.utils import parse_phonopy_force_constants
from hilde.harmonic_analysis.dynamical_matrix import get_frequencies

parent = Path(__file__).parent
assets = parent / "assets_remap"
uc_filename = assets / "geometry.in.primitive"
sc_filename = assets / "geometry.in.supercell"
fc_filename = assets / "FORCE_CONSTANTS"

frequencies = np.loadtxt(assets / "frequencies.dat")


def test_remap():

    atoms = read(sc_filename, format="aims")

    fc = parse_phonopy_force_constants(
        uc_filename=uc_filename,
        sc_filename=sc_filename,
        fc_filename=fc_filename,
        two_dim=True,
    )

    freqs = get_frequencies(fc, masses=atoms.get_masses())

    assert np.linalg.norm(freqs - frequencies) < 1e-12


if __name__ == "__main__":
    test_remap()
