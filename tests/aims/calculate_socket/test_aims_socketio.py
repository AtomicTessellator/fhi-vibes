"""test the MD workflow"""
import subprocess as sp
from pathlib import Path

import numpy as np
import pytest

from vibes.trajectory import reader


parent = Path(__file__).parent

run_command = "vibes run singlepoint aims.in"

task_input_file = parent / "aims.in"
geometry_input_files = sorted(parent.glob("geometry.in.???"))

trajectory_ref_file = parent / "trajectory.ref.son"


@pytest.mark.parametrize("socketio", (True, False))
def test_aims_calculate(tmp_path, socketio):
    # read input file
    task_input = open(task_input_file).read().replace("SOCKETIO", str(socketio))
    (tmp_path / task_input_file.name).write_text(task_input)

    for file in geometry_input_files:
        (tmp_path / file.name).symlink_to(file)

    sp.run(run_command.split(), cwd=tmp_path, check=True)

    old_traj = reader(trajectory_ref_file)
    new_traj = reader(tmp_path / "aims" / "trajectory.son")

    _check_ref(old_traj, new_traj)


def _check_ref(trajectory1, trajectory2):
    for (a0, a1) in zip(trajectory1, trajectory2):
        r0, r1 = a0.calc.results, a1.calc.results
        for key in r0:
            if key in r1:
                np.testing.assert_allclose(r0[key], r1[key], rtol=1e-3)
