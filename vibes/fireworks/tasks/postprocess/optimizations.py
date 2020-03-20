"""Postprocess steps for k-grid optimizations"""
from pathlib import Path

from vibes.helpers.fileformats import last_from_yaml


def load_last_step(atoms, calculator, workdir, trajectory_file):
    """Loads the last step from a trajectory and update returns calculator objects"""
    trajectory_file = Path(workdir) / trajectory_file
    last_step_dict = last_from_yaml(trajectory_file)

    for key, val in last_step_dict["atoms"].items():
        atoms[key] = val

    calculator["results"] = last_step_dict["calculator"]

    return trajectory_file, atoms, calculator


def move_trajectory_file(trajectory_file):
    """Move a trajectory to a new file name"""
    split_trajectory_file = trajectory_file.split(".")
    try:
        temp_list = split_trajectory_file[-2].split("_")
        temp_list[-1] = str(int(temp_list[-1]) + 1)
        split_trajectory_file[-2] = "_".join(temp_list)
        trajectory_file = ".".join(split_trajectory_file)
    except ValueError:
        split_trajectory_file[-2] += "_restart_1"
        trajectory_file = ".".join(split_trajectory_file)
