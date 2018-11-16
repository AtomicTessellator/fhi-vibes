""" read YAML trajectories """

import json
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from hilde.helpers.fileformats import from_yaml
from hilde.helpers.converters import dict2results


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

    trajectory = []
    for obj in pre_trajectory:
        # Atoms
        atoms_dict = obj['atoms']
        calc_dict = {**calculator_data, **obj['calculator']}

        atoms = dict2results(atoms_dict, calc_dict)

        trajectory.append(atoms)

    return trajectory
