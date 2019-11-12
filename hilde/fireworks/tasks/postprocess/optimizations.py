"""Postprocess steps for k-grid optimizations"""
from pathlib import Path
from hilde.helpers.fileformats import last_from_yaml


def load_last_step(atoms, calc, workdir, trajectory):
    """Loads the last step from a trajectory and update the atoms and calculator objects"""
    trajectory = Path(workdir) / trajectory
    last_step_dict = last_from_yaml(trajectory)

    for key, val in last_step_dict["atoms"].items():
        atoms[key] = val

    calc["results"] = last_step_dict["calculator"]

    return trajectory, atoms, calc


def move_trajectory_file(trajectory):
    """Move a trajectory to a new file name"""
    new_traj_list = trajectory.split(".")
    try:
        temp_list = new_traj_list[-2].split("_")
        temp_list[-1] = str(int(temp_list[-1]) + 1)
        new_traj_list[-2] = "_".join(temp_list)
        trajectory = ".".join(new_traj_list)
    except ValueError:
        new_traj_list[-2] += "_restart_1"
        trajectory = ".".join(new_traj_list)
