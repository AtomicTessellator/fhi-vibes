#!/usr/bin/env python
# coding: utf-8

from pathlib import Path

from vibes.trajectory import Trajectory, analysis

parent = Path(__file__).parent

file = parent / "mgo.son"


def test_summary(file=file):
    trajectory = Trajectory.read(file)
    analysis.summary(trajectory.dataset, plot=True)
    assert Path("md_summary.pdf").exists()


if __name__ == "__main__":
    test_summary()
