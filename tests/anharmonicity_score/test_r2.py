"""test the anharmonicity quantification"""

from pathlib import Path
import numpy as np

from hilde.trajectory import reader
from hilde.tdep.wrapper import parse_tdep_remapped_forceconstant
import hilde.anharmonicity_score as score


parent = Path(__file__).parent


def test_r2():
    fc = parse_tdep_remapped_forceconstant(parent / "outfile.forceconstant_remapped")

    trajectory = reader(parent / "trajectory.son")
    trajectory.force_constants = fc

    supercell = trajectory.supercell

    f_dft, f_ha = trajectory.forces, trajectory.forces_harmonic

    r2 = score.get_r2(f_dft, f_ha)
    ref_r2 = 0.56517063272
    match_r2 = np.allclose(r2, ref_r2, rtol=1e-8)

    assert match_r2, (r2, ref_r2)

    r2_per_atom = score.get_r2_per_atom(f_dft, f_ha, supercell, reduce_by_symmetry=True)
    ref_r2_per_atom = [0.316842934892, 0.544809575396, 0.580653088592]
    match_r2_per_atom = np.allclose(r2_per_atom, ref_r2_per_atom, rtol=1e-8)

    assert match_r2_per_atom, (r2_per_atom, ref_r2_per_atom)


if __name__ == "__main__":
    test_r2()
