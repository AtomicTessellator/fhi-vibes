"""test the anharmonicity quantification"""

from pathlib import Path
import numpy as np

from vibes.trajectory import reader
from vibes.tdep.wrapper import parse_tdep_remapped_forceconstant
import vibes.anharmonicity_score as score


parent = Path(__file__).parent

fc = parse_tdep_remapped_forceconstant(parent / "outfile.forceconstant_remapped")

trajectory = reader(parent / "trajectory.son")
trajectory.set_force_constants_remapped(fc)


# def test_r2():
#     supercell = trajectory.supercell
#
#     f_dft, f_ha = trajectory.forces, trajectory.forces_harmonic
#
#     r2 = score.get_r2(f_dft, f_ha)
#     ref_r2 = 0.56517063272
#     match_r2 = np.allclose(r2, ref_r2, rtol=1e-8)
#
#     assert match_r2, (r2, ref_r2)
#
#     r2_per_atom = score.get_r2_per_atom(f_dft, f_ha, supercell, reduce_by_symmetry=True)
#     ref_r2_per_atom = [0.316842934892, 0.544809575396, 0.580653088592]
#     match_r2_per_atom = np.allclose(r2_per_atom, ref_r2_per_atom, rtol=1e-8)
#
#     assert match_r2_per_atom, (r2_per_atom, ref_r2_per_atom)


def test_sigma():
    df = score.get_dataframe(trajectory.dataset)

    assert np.allclose(df.sigma, 0.659416), df.sigma

    sigmas = (df["sigma [Cs]"], df["sigma [Pb]"], df["sigma [I]"])
    sigma_per_atom = [float(v) for v in sigmas]
    ref_sigma_per_atom = [0.826162, 0.674311, 0.647513]
    match_sigma_per_atom = np.allclose(sigma_per_atom, ref_sigma_per_atom)

    assert match_sigma_per_atom, (sigma_per_atom, ref_sigma_per_atom)


def test_sigma_mode():
    series = score.get_sigma_per_mode(trajectory.dataset)

    assert np.allclose(series.mean(), 1.099467483), series.mean()


if __name__ == "__main__":
    test_sigma()
    test_sigma_mode()
