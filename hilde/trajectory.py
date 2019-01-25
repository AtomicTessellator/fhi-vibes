""" tools for storing MD trajectories

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import json

from hilde import __version__ as version
from hilde.helpers.converters import results2dict, dict2results
from hilde.helpers.fileformats import to_yaml, from_yaml
from hilde.helpers.hash import hash_atoms


def step2file(atoms, calc, file="trajectory.yaml", append_cell=False):
    """ Save the current step """

    to_yaml(results2dict(atoms, calc, append_cell), file)


def metadata2file(metadata, file="metadata.yaml"):
    """ save metadata to file """

    if metadata is None:
        metadata = {}

    to_yaml({**metadata, "hilde": {"version": version}}, file, mode="w")


def get_hashes_from_trajectory(trajectory):
    """ return all hashes from trajectory """

    try:
        traj = reader(trajectory)
    except (FileNotFoundError, KeyError):
        return []

    hashes = []
    for atoms in traj:
        try:
            hashes.append(atoms.info["hash"])
        except (KeyError, AttributeError):
            hashes.append(hash_atoms(atoms))

    return hashes


def reader(file, get_metadata=False):
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

    if "MD" in metadata:
        md_metadata = metadata["MD"]

    trajectory = []
    for obj in pre_trajectory:

        atoms_dict = {**pre_atoms_dict, **obj["atoms"]}

        # remember that the results need to go to a dedicated results dict in calc
        calc_dict = {**pre_calc_dict, "results": obj["calculator"]}

        atoms = dict2results(atoms_dict, calc_dict)

        # info
        if "MD" in metadata:
            if "dt" in atoms.info:
                atoms.info["dt"] /= md_metadata["fs"]
        elif "info" in obj:
            info = obj["info"]
            atoms.info.update(info)

        trajectory.append(atoms)
    if get_metadata:
        return trajectory, metadata
    return trajectory
