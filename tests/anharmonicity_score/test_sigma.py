"""test the anharmonicity quantification"""

from pathlib import Path

import numpy as np

import vibes.anharmonicity_score as score
from vibes.trajectory import reader

parent = Path(__file__).parent

trajectory = reader(parent / "trajectory.nc")


def test_sigma():
    df = score.get_dataframe(trajectory.dataset)

    assert np.allclose(df.sigma, 0.659416), df.sigma

    sigmas = (df["sigma [Cs]"], df["sigma [Pb]"], df["sigma [I]"])
    sigma_per_atom = [float(v) for v in sigmas]
    ref_sigma_per_atom = [0.826162, 0.674311, 0.647513]
    match_sigma_per_atom = np.allclose(sigma_per_atom, ref_sigma_per_atom)

    assert match_sigma_per_atom, (sigma_per_atom, ref_sigma_per_atom)


def test_sigma_mode():
    _ = score.get_sigma_per_mode(trajectory.dataset)

    # numerically not stable!!!
    # assert np.allclose(series.iloc[3:].mean(), 1.099467483), series.mean()


if __name__ == "__main__":
    test_sigma()
    test_sigma_mode()
