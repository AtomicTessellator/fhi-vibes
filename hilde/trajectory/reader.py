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

    pre_calc_dict = metadata["calculator"]
    pre_atoms_dict = metadata["atoms"]

    if "numbers" in pre_atoms_dict and "symbols" in pre_atoms_dict:
        del pre_atoms_dict["symbols"]

    trajectory = []
    for obj in pre_trajectory:

        atoms_dict = {**pre_atoms_dict, **obj["atoms"]}

        # remember that the results need to go to a dedicated results dict in calc
        calc_dict = {**pre_calc_dict, "results": obj["calculator"]}

        atoms = dict2results(atoms_dict, calc_dict)

        trajectory.append(atoms)

    return trajectory
