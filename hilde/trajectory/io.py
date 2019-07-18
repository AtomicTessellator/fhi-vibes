"""Trajectory File I/O"""

import json
from hilde import son
from hilde.helpers.converters import dict2atoms
from hilde.helpers import Timer, warn, talk
from hilde.trajectory import Trajectory
from hilde.helpers.utils import progressbar


def reader(file="trajectory.son", get_metadata=False, verbose=True):
    """Convert information in trajectory and metadata files to atoms objects and return them

    Parameters
    ----------
    trajectory: str
        Trajectory file to pull the structures from
    get_metadata: bool
        If True return the metadata
    verbose: bool
        If True print more information to the screen

    Returns
    -------
    trajectory: Trajectory
        The trajectory from the file
    metadata: dict
        The metadata for the trajectory
    """

    timer = Timer(f"Parse trajectory")

    try:
        metadata, pre_trajectory = son.load(file, verbose=verbose)
    except json.decoder.JSONDecodeError:
        metadata, pre_trajectory = son.load(file, verbose=verbose)

    # legacy of trajectory.yaml
    if metadata is None:
        msg = f"metadata in {file} appears to be empty, assume old convention w/o === "
        msg += f"was used. Let's see"
        warn(msg, level=1)
        metadata = pre_trajectory.pop(0)

    pre_calc_dict = metadata["calculator"]
    pre_atoms_dict = metadata["atoms"]

    if "numbers" in pre_atoms_dict and "symbols" in pre_atoms_dict:
        del pre_atoms_dict["symbols"]

    if "MD" in metadata:
        md_metadata = metadata["MD"]

    if not pre_trajectory:
        if get_metadata:
            talk(".. trajectory empty, return ([], metadata)")
            return [], metadata
        talk(".. trajectory empty, return []")
        return []

    trajectory = Trajectory(metadata=metadata)
    prefix = ".. create atoms: "
    for obj in progressbar(pre_trajectory, prefix=prefix):

        atoms_dict = {**pre_atoms_dict, **obj["atoms"]}

        # remember that the results need to go to a dedicated results dict in calc
        calc_dict = {**pre_calc_dict, "results": obj["calculator"]}

        atoms = dict2atoms(atoms_dict, calc_dict)

        # info
        if "MD" in metadata:
            if "dt" in atoms.info:
                atoms.info["dt_fs"] = atoms.info["dt"] / md_metadata["fs"]
        elif "info" in obj:
            info = obj["info"]
            atoms.info.update(info)

        # compatibility with older trajectories
        if "MD" in obj:
            atoms.info.update(obj["MD"])

        # preserve metadata
        if "metadata" in obj:
            atoms.info.update({"metadata": obj["metadata"]})

        trajectory.append(atoms)

    timer("done")

    if get_metadata:
        return trajectory, metadata
    return trajectory
