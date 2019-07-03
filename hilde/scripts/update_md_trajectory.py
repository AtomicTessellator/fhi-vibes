""" Update trajectory files of old format """

from argparse import ArgumentParser
import shutil

from hilde.io import read
from hilde.trajectory import reader
from hilde.helpers import talk


def update_trajectory(trajectory, uc=None, sc=None, format="aims"):
    """update TRAJECTORY by adding unit cell and supercell"""
    traj = reader(trajectory)
    new_trajectory = "temp.son"

    if uc:
        atoms = read(uc, format=format)
        traj.primitive = atoms

    if sc:
        atoms = read(sc, format=format)
        traj.supercell = atoms

    traj.write(file=new_trajectory)

    fname = f"{trajectory}.bak"
    talk(f".. back up old trajectory to {fname}")
    shutil.copy(trajectory, fname)
    talk(f".. write new trajectory to {trajectory}")
    shutil.move(new_trajectory, trajectory)


def main():
    """ main routine """
    parser = ArgumentParser(description="Update trajectory file")
    parser.add_argument("trajectory")
    parser.add_argument("-uc", help="Add a (primitive) unit cell")
    parser.add_argument("-sc", help="Add the respective supercell")
    parser.add_argument("--format", default="aims")
    args = parser.parse_args()

    update_trajectory(args.trajectory, args.uc, args.sc, args.format)


if __name__ == "__main__":
    main()
