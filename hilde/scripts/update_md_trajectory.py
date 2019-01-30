""" Update trajectory files of old format """

from argparse import ArgumentParser
import shutil

from hilde.helpers.fileformats import from_yaml, to_yaml


def main():
    """ main routine """
    parser = ArgumentParser(description="Update trajectory file")
    parser.add_argument("trajectory")
    args = parser.parse_args()

    trajectory = args.trajectory
    new_trajectory = "temp.yaml"

    metadata, *traj = from_yaml(trajectory)

    to_yaml(metadata, new_trajectory, mode="w")

    for elem in traj:
        if "MD" in elem:
            info = elem.pop("MD")
            elem["atoms"].update({"info": info})
        to_yaml(elem, new_trajectory)

    shutil.copy(trajectory, f"{trajectory}.bak")
    shutil.move(new_trajectory, trajectory)


if __name__ == "__main__":
    main()

