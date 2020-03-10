""" Compute the phonon fingerprints for supercells of different size """

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.dft.kpoints import get_special_points

from vibes.helpers.supercell import make_cubic_supercell
from vibes.materials_fp.material_fingerprint import (
    get_phonon_bs_fp,
    get_phonon_dos_fp,
    scalar_product,
)
from vibes.phonopy import wrapper as ph


atoms = bulk("Al")

# Space group information
special_points = get_special_points(atoms.cell)

# Calculator setup

# conventional supercell matrix
cmatrix = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])


def test_all():
    # run phonon calculation for several supercell sizes and compute fingerprints
    fps_bs = []
    fps_dos = []
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
        phonon.run_mesh([15, 15, 15])
        phonon.run_total_dos(freq_min=0.0, freq_max=10.0, freq_pitch=0.01)

        fp_bs = get_phonon_bs_fp(phonon, special_points, binning=False)[0]
        fps_bs.append(fp_bs)

        fp_dos = get_phonon_dos_fp(phonon, nbins=1001)
        fps_dos.append(fp_dos)

    fps_bs = np.asarray(fps_bs)

    fp_diffs = abs(fps_bs - fps_bs[-1]).max(axis=2)

    sps = [
        0.0,
        scalar_product(
            fps_dos[1], fps_dos[0], col=1, pt=0, normalize=False, tanimoto=True
        ),
        scalar_product(
            fps_dos[2], fps_dos[1], col=1, pt=0, normalize=False, tanimoto=True
        ),
    ]
    print(
        "n_atoms   "
        + " ".join([f"{k:9s}" for k in special_points.keys()])
        + "similarity_score"
    )
    for nn, fp, sp in zip(n_atoms, fp_diffs, sps):
        print(f"{nn:4d}: " + " ".join([f"{f:9.3e}" for f in fp]) + f"    {sp:9.3e}")

    assert all(3 < fp < 9 for fp in fps_bs[-1][1])


if __name__ == "__main__":
    test_all()
