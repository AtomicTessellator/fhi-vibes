"""
Script to compute supercells from inputself.
Similar to generate_structure from TDEP.
"""
from argparse import ArgumentParser as argpars
import numpy as np
from hilde.parsers import read_structure
from hilde.helpers.supercell import make_cubic_supercell, make_supercell
from hilde.helpers.geometry import get_cubicness


def print_matrix(matrix):
    rep = [" [{}, {}, {}]".format(*elem) for elem in matrix]
    # join to have a comma separated list
    rep = f",\n".join(rep)
    # add leading [ and trailing ]
    rep = f"[{rep[1:]}"
    rep += "]"
    print(rep)


def main():
    parser = argpars(description="Read geometry create supercell")
    parser.add_argument("geom", type=str, help="geometry input file")
    parser.add_argument("-n", type=int, help="target size")
    parser.add_argument("-d", type=int, nargs=9, help="supercell matrix")
    parser.add_argument("--dry", action="store_true", help="Do not write output file")
    parser.add_argument("--format", default="aims")
    args = parser.parse_args()

    fname = args.geom
    print(f"Find supercell for")
    cell = read_structure(fname, format=args.format)
    cell.inform(fname=fname)

    print("\nSettings:")
    if args.n:
        print(f"  Target number of atoms: {args.n}")
        supercell, smatrix = make_cubic_supercell(cell, args.n)
        print(f"  Found number of atoms:  {len(supercell)}")
    elif args.d:
        smatrix = np.array(args.d).reshape(3, 3)
        supercell = make_supercell(cell, smatrix)
    else:
        exit("Please specify either a target cell size or a supercell matrix")

    print(f"\nSupercell matrix:")
    print_matrix(smatrix)
    print(f"\nSuperlattice:")
    print(supercell.cell)
    print(
        f"Cubicness:  {get_cubicness(supercell.cell):.3f} "
        + f"({get_cubicness(supercell.cell)**3:.3f})"
    )

    if not args.dry:
        output_filename = f"{args.geom}.supercell"
        supercell.write(output_filename, scaled=False)
        print(f"Supercell written to {output_filename}")


if __name__ == "__main__":
    main()
