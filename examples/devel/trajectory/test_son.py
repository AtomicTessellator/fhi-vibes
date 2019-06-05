#!/usr/bin/env python
# coding: utf-8

import son
from hilde.helpers.fileformats import dict2json

new_trajectory = "test.son"

meta, traj = son.load("trajectory.son")

son.dump(meta, new_trajectory, dumper=dict2json, is_metadata=True)

for atoms in traj:
    print(atoms)
    son.dump(atoms, new_trajectory, dumper=dict2json)


new_meta, new_traj = son.load(new_trajectory)

assert open("trajectory.son").read() == open("test.son").read()


assert abs(new_meta["MD"]["timestep"] - meta["MD"]["timestep"]) < 1e-14
