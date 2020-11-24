"""test the MD workflow"""
import subprocess as sp
from pathlib import Path

import numpy as np
import pytest

from vibes.trajectory import reader


parent = Path(__file__).parent

run_command = "vibes run md"

md_input_files = (parent / file for file in ("md.in", "md.in.nosocket"))
geometry_input_file = parent / "geometry.in"

trajectory_ref_file = parent / "trajectory.ref.son"


@pytest.mark.parametrize("md_input_file", md_input_files)
def test_aims_md(tmp_path, md_input_file):
    (tmp_path / "geometry.in").symlink_to(geometry_input_file)
    (tmp_path / "md.in").symlink_to(md_input_file)

    sp.run(run_command.split(), cwd=tmp_path)

    old_traj = reader(trajectory_ref_file)
    new_traj = reader(tmp_path / "md" / "trajectory.son")

    _check_ref(old_traj, new_traj)


def _check_ref(trajectory1, trajectory2):
    for (a0, a1) in zip(trajectory1, trajectory2):
        r0, r1 = a0.calc.results, a1.calc.results
        for key in r0:
            if key in r1:
                np.testing.assert_allclose(r0[key], r1[key], rtol=1e-3)
