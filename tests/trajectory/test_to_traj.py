#!/usr/bin/env python
# coding: utf-8

from pathlib import Path

from ase.io import Trajectory

from vibes.helpers.hash import hash_atoms, hash_atoms_and_calc
from vibes.trajectory import reader

parent = Path(__file__).parent
traj = reader(parent / "trajectory.son", single_point_calculator=False)

traj_file = Path(parent / "test.traj")


def test_to_traj():

    traj.to_traj(traj_file)

    ase_traj = Trajectory(traj_file)

    for i, atoms in enumerate(traj):
        atoms_from_ase = ase_traj[i]

        assert hash_atoms(atoms, velocities=True) == hash_atoms(
            atoms_from_ase, velocities=True
        )

        assert hash_atoms_and_calc(atoms, ignore_results=False) == hash_atoms_and_calc(
            atoms_from_ase, ignore_results=False
        )

    traj_file.unlink()


if __name__ == "__main__":
    test_to_traj()
