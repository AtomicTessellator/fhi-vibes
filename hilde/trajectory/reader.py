""" read YAML trajectories """

from pathlib import Path
from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from hilde.helpers.fileformats import from_json
from hilde.helpers.fileformats import from_yaml


def reader(file, metadata):
    """ convert information in trajectory and metadata files to atoms objects
     and return them """

    if str(file).endswith("json"):
        pre_trajectory = from_json(file)
    elif str(file).endswith("yaml"):
        pre_trajectory = from_yaml(file)
    else:
        raise Exception("Only json and yaml supported")

    if str(metadata).endswith("json"):
        metadata = from_json(metadata)
    elif str(metadata).endswith("yaml"):
        metadata = from_yaml(metadata)
    else:
        raise Exception("Only json and yaml supported")

    calculator_data = metadata["calculator"]

    pre_atoms_dict = metadata["atoms"]

    # only one of symbols or numbers is needed
    del pre_atoms_dict["symbols"]

    pbc = False
    if "cell" in pre_atoms_dict:
        pbc = True

    trajectory = []
    for obj in pre_trajectory:

        velocities = obj["atoms"].pop("velocities")
        atoms_dict = {**pre_atoms_dict, **obj["atoms"]}
        atoms = Atoms(**atoms_dict, pbc=pbc)
        atoms.set_velocities(velocities)

        results = obj["calculator"]
        calc = SinglePointCalculator(atoms, **results)
        calc.name = calculator_data["name"]
        calc.parameters.update(calculator_data["params"])

        atoms.calc = calc

        trajectory.append(atoms)

    return trajectory
