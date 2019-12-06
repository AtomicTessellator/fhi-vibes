#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
import numpy as np
from hilde.trajectory import Trajectory

parent = Path(__file__).parent

traj0 = Trajectory.read(parent / "trajectory.son")
traj0.write(parent / "trajectory.nc")
traj1 = Trajectory.read(parent / "trajectory.nc")


def test_compare_son_nc(traj0=traj0, traj1=traj1):
    for a0, a1 in zip(traj0, traj1):
        assert np.allclose(a0.positions, a1.positions)
        assert np.allclose(a0.get_velocities(), a1.get_velocities())
        assert np.allclose(a0.get_forces(), a1.get_forces())
        assert np.allclose(a0.get_total_energy(), a1.get_total_energy())
        assert np.allclose(a0.get_kinetic_energy(), a1.get_kinetic_energy())


if __name__ == "__main__":
    test_compare_son_nc()
