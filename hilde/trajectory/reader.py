""" read YAML trajectories """
import json
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from hilde.helpers.fileformats import from_yaml
from hilde.helpers.converters import dict2results


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
        info = {}
        if "MD" in obj:
            info = obj["MD"]
            info["dt"] /= md_metadata["fs"]
        elif "info" in obj:
            info = obj["info"]

        atoms.info = info

        trajectory.append(atoms)
    if get_metadata:
        return trajectory, CaseInsensitiveDict(metadata)
    return trajectory

class CaseInsensitiveDict(dict):
    def __init__(self, d):
        self._d = d
        self._s = dict((k.lower(), k) for k in d)
    def __contains__(self, k):
        return k.lower() in self._s
    def __len__(self):
        return len(self._s)
    def __iter__(self):
        return iter(self._s)
    def __getitem__(self, k):
        return self._d[self._s[k.lower()]]
