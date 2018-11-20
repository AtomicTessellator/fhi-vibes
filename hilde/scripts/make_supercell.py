"""
Script to compute supercells from inputself.
Similar to generate_structure from TDEP.
"""
from argparse import ArgumentParser as argpars
import numpy as np
from ase.io import read
from hilde.structure.io import inform
from hilde.helpers.supercell import make_cubic_supercell, make_supercell
from hilde.helpers.geometry import get_cubicness
from hilde.helpers.maths import get_3x3_matrix


def print_matrix(matrix, indent=2):
    ind = indent * " "
    rep = [" [{}, {}, {}]".format(*elem) for elem in matrix]
    # join to have a comma separated list
    rep = f",\n{ind}".join(rep)
    # add leading [ and trailing ]
    rep = f"{ind}[{rep[1:]}"
    rep += "]"
    print(rep)


def main():
    parser = argpars(description="Read geometry create supercell")
    parser.add_argument("geom", type=str, help="geometry input file")
    parser.add_argument("-n", type=int, help="target size")
    parser.add_argument("-d", type=int, nargs="+", help="supercell matrix")
    parser.add_argument("--deviation", type=float, default=0.2)
    parser.add_argument("--dry", action="store_true", help="Do not write output file")
    parser.add_argument("--format", default="aims")
    args = parser.parse_args()

    fname = args.geom
    print(f"Find supercell for")
    cell = read(fname, format=args.format)
    inform(cell)

    if args.n:
        print("\nSettings:")
        print(f"  Target number of atoms: {args.n}")
        supercell, smatrix = make_cubic_supercell(
            cell, args.n, deviation=args.deviation
        )
    elif args.d:
        smatrix = get_3x3_matrix(args.d)
        supercell = make_supercell(cell, smatrix)
    else:
        exit("Please specify either a target cell size or a supercell matrix")

    print(f"\nSupercell matrix:")
    print(" 1d: {}".format(np.array2string(smatrix.flatten(), separator=", ")))
    print(" 2d:")
    print_matrix(smatrix, indent=0)

    print(f"\nNumber of atoms:  {len(supercell)}")

    print(f"\nSuperlattice:")
    print(supercell.cell)
    print(
        f"Cubicness:  {get_cubicness(supercell.cell):.3f} "
        + f"({get_cubicness(supercell.cell)**3:.3f})"
    )

    if not args.dry:
        output_filename = f"{args.geom}.supercell"
        supercell.write(output_filename, scaled=False, format=args.format)
        print(f"Supercell written to {output_filename}")


if __name__ == "__main__":
    main()
