"""test the anharmonicity quantification"""

import numpy as np

from hilde.trajectory import reader
from hilde.tdep.wrapper import parse_tdep_forceconstant

import hilde.anharmonicity_score as score


def test_r2():
    fc = parse_tdep_forceconstant("outfile.forceconstant_remapped")

    trajectory = reader("trajectory.son")

    supercell = trajectory.supercell

    f_dft, f_ha = score.get_forces_from_trajectory(trajectory, supercell, fc)

    r2 = score.get_r2(f_dft, f_ha)

    _, r2_per_atom = score.get_r2_per_atom(
        f_dft, f_ha, supercell, reduce_by_symmetry=True
    )

    ref_r2 = 0.56517063272
    ref_r2_per_atom = [0.316842934892, 0.544809575396, 0.580653088592]

    match_r2 = np.allclose(r2, ref_r2, rtol=1e-8)
    match_r2_per_atom = np.allclose(r2_per_atom, ref_r2_per_atom, rtol=1e-8)

    assert match_r2, r2
    assert match_r2_per_atom, r2_per_atom


if __name__ == "__main__":
    test_r2()