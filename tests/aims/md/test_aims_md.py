"""test the MD workflow"""
import subprocess as sp
from pathlib import Path

import numpy as np
import pytest

from vibes.trajectory import reader
from vibes.helpers.stresses import has_stresses


parent = Path(__file__).parent

run_command = "vibes run md"

md_input_files = (parent / file for file in ("md.in", "md.in.nosocket"))
geometry_input_file = parent / "geometry.in"

trajectory_ref_file = parent / "trajectory.ref.son"


@pytest.mark.parametrize("md_input_file", md_input_files)
def test_aims_md(tmp_path, md_input_file):
    (tmp_path / "geometry.in").symlink_to(geometry_input_file)

    md_in = tmp_path / "md.in"

    md_in.symlink_to(md_input_file)
    sp.run(run_command.split(), cwd=tmp_path)
    md_in.unlink()
    # test restart
    md_in.symlink_to(md_input_file.with_suffix(md_input_file.suffix + ".2"))
    sp.run(run_command.split(), cwd=tmp_path)

    old_traj = reader(trajectory_ref_file)
    new_traj = reader(tmp_path / "md" / "trajectory.son")

    _check_ref(old_traj, new_traj)


def _check_ref(trajectory1, trajectory2):
    for i in range(len(trajectory1)):
        # we can't use zip() because it truncates if the trajectories
        # have different lengths
        a0, a1 = trajectory1[i], trajectory2[i]
        r0, r1 = a0.calc.results, a1.calc.results
        assert has_stresses(a0) == has_stresses(a1)
        for key in r0:
            if key in r1:
                np.testing.assert_allclose(r0[key], r1[key], rtol=1e-3)
