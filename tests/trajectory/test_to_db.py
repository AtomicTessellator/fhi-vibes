#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
from vibes.trajectory import reader
from vibes.helpers.hash import hash_atoms, hash_atoms_and_calc
from ase.db import connect

parent = Path(__file__).parent
traj = reader(parent / "trajectory.son", single_point_calc=False)

db_file = Path(parent / "db.json")

traj.to_db(db_file)

db = connect(db_file)

for i, atoms in enumerate(traj):
    atoms_from_db = db.get_atoms(
        i+1, attach_calculator=True, add_additional_information=True
    )

    assert hash_atoms(atoms, velocities=True) == hash_atoms(
        atoms_from_db, velocities=True
    )

    assert hash_atoms_and_calc(atoms, ignore_results=False) == hash_atoms_and_calc(
        atoms_from_db, ignore_results=False
    )
    assert atoms.info == atoms_from_db.info["data"]["info"]

db_file.unlink()
