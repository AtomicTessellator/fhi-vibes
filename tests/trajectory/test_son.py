#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
from hilde import son

parent = Path(__file__).parent
new_trajectory = parent / "test.son"

meta, traj = son.load(parent / "trajectory.son")

son.dump(meta, new_trajectory, is_metadata=True)

for atoms in traj:
    son.dump(atoms, new_trajectory)


new_meta, new_traj = son.load(new_trajectory)

assert open(parent / "trajectory.son").read() == open(new_trajectory).read()

assert abs(new_meta["MD"]["timestep"] - meta["MD"]["timestep"]) < 1e-14

# clean up
new_trajectory.unlink()
