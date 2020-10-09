#!/usr/bin/env python
# coding: utf-8

from pathlib import Path

from vibes import keys
from vibes.trajectory import Trajectory


parent = Path(__file__).parent

file = parent / "mgo.son"


def test_parse_trajectory(file=file):
    trajectory = Trajectory.read(file)
    trajectory.compute_heat_fluxes_from_stresses()
    DS = trajectory.dataset

    for key in (
        keys.heat_flux,
        keys.heat_flux_aux,
    ):
        assert key in DS


if __name__ == "__main__":
    test_parse_trajectory()
