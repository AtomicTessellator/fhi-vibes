""" Compute the phonon fingerprints for supercells of different size """

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.dft.kpoints import get_special_points

from vibes.helpers.supercell import make_cubic_supercell
from vibes.materials_fp.material_fingerprint import get_phonon_bs_fingerprint_phononpy
from vibes.phonopy import wrapper as ph

atoms = bulk("Al")

# Space group information
special_points = get_special_points(atoms.cell)

# Calculator setup

# conventional supercell matrix
cmatrix = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])


def test_all():
    # run phonon calculation for several supercell sizes and compute fingerprints
    fps = []
    n_atoms = []
    for nn in [4, 32, 108]:

        # smatrix = a * cmatrix
        supercell, smatrix = make_cubic_supercell(atoms, nn)

        n_a = len(supercell)
        print(f"compute for {n_a} atoms")
        n_atoms.append(n_a)

        phonon, sc, scs = ph.preprocess(atoms, smatrix.T)

        force_sets = []
        for cell in scs:
            cell.calc = EMT()
            force_sets.append(cell.get_forces())

        phonon.produce_force_constants(force_sets)

        # REM: binning=False is optional
        fp = get_phonon_bs_fingerprint_phononpy(phonon, special_points, binning=False)[
            0
        ]
        fps.append(fp)
    fps = np.asarray(fps)

    # Compute difference to largest supercell and choose largest deviation at each
    # q point
    fp_diffs = abs(fps - fps[-1]).max(axis=2)

    print("n_atoms   " + " ".join([f"{k:9s}" for k in special_points.keys()]))
    for nn, fp in zip(n_atoms, fp_diffs):
        print(f"{nn:4d}: " + " ".join([f"{f:9.3e}" for f in fp]))

    assert all(3 < fp < 9 for fp in fps[-1][1])


if __name__ == "__main__":
    test_all()
