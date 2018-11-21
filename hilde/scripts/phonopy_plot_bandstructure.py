""" Take a pickled phonopy object and print bandstructure """

from argparse import ArgumentParser
from hilde.helpers.pickletools import pread
from hilde.phonopy.wrapper import plot_bandstructure


def main():
    """ main routine """
    parser = ArgumentParser(description="plot band structure from a phonopy object")
    parser.add_argument("file", help="pickled phonopy object")
    parser.add_argument("--outfile", default="bandstructure.pdf")
    args = parser.parse_args()

    phonon = pread(args.file)

    plot_bandstructure(phonon, file=args.outfile)


if __name__ == "__main__":
    main()
