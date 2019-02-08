""" Update trajectory files of old format """

from argparse import ArgumentParser

from hilde.trajectory import reader


def main():
    """ main routine """
    parser = ArgumentParser(description="Update trajectory file")
    parser.add_argument("trajectory")
    parser.add_argument("-s", "--skip", default=1, type=int)
    parser.add_argument("-f", "--folder", default=".")
    parser.add_argument("-uc", '--ucposcar', action='store_true')
    args = parser.parse_args()

    trajectory = reader(args.trajectory)

    trajectory.to_tdep(folder=args.folder, skip=args.skip, ucposcar=args.ucposcar)


if __name__ == "__main__":
    main()
