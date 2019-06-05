#!/usr/bin/env python
# coding: utf-8

from hilde import son

new_trajectory = "test.son"

meta, traj = son.load("trajectory.son")

son.dump(meta, new_trajectory, is_metadata=True)

for atoms in traj:
    son.dump(atoms, new_trajectory)


new_meta, new_traj = son.load(new_trajectory)

assert open("trajectory.son").read() == open("test.son").read()


assert abs(new_meta["MD"]["timestep"] - meta["MD"]["timestep"]) < 1e-14
