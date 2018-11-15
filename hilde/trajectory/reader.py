""" read YAML trajectories """

import json
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from hilde.helpers.fileformats import from_yaml


def reader(file):
    """ convert information in trajectory and metadata files to atoms objects
     and return them """

    try:
        metadata, *pre_trajectory = from_yaml(file, use_json=True)
    except json.decoder.JSONDecodeError:
        metadata, *pre_trajectory = from_yaml(file, use_json=False)

    calculator_data = metadata["calculator"]
    pre_atoms_dict = metadata["atoms"]

    if "numbers" in pre_atoms_dict and "symbols" in pre_atoms_dict:
        del pre_atoms_dict["symbols"]

    pbc = False
    if "cell" in pre_atoms_dict:
        pbc = True

    trajectory = []
    for obj in pre_trajectory:
        # Atoms
        try:
            velocities = obj["atoms"].pop("velocities")
        except KeyError:
            velocities = None

        atoms_dict = {**pre_atoms_dict, **obj["atoms"]}
        atoms = Atoms(**atoms_dict, pbc=pbc)

        if velocities is not None:
            atoms.set_velocities(velocities)

        # Calculator
        results = obj["calculator"]
        calc = SinglePointCalculator(atoms, **results)
        calc.name = calculator_data["name"]
        calc.parameters.update(calculator_data["params"])

        atoms.calc = calc

        trajectory.append(atoms)

    return trajectory
