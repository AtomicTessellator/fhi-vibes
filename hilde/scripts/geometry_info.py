import os
from argparse import ArgumentParser as argpars
from hilde.parsers import read_structure


def main():
    """ print geometry information """
    parser = argpars(description="Read geometry and print some symmetry info")
    parser.add_argument("geom", type=str, help="geometry input file")
    parser.add_argument(
        "-t",
        "--tolerance",
        type=float,
        default=1e-5,
        help="symmetry tolerance (symprec)",
    )
    parser.add_argument(
        "--format", default="aims", help="which format should ASE use to read file?"
    )
    args = parser.parse_args()

    ### Greet
    fname = args.geom
    cell = read_structure(fname, symprec=args.tolerance, format=args.format)
    cell.inform(fname=fname)


if __name__ == "__main__":
    qlaunch()
