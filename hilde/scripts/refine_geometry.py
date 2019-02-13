from argparse import ArgumentParser as argpars

from hilde.io import read, write, inform
from hilde.spglib.wrapper import refine_cell, standardize_cell


def main():
    parser = argpars(description="Read geometry and use spglib to refine")
    parser.add_argument("geometry")
    parser.add_argument("--prim", action="store_true", help="store primitive cell")
    parser.add_argument("--conv", action="store_true", help="store conventional cell")
    parser.add_argument("-t", "--tolerance", type=float, default=1e-5)
    parser.add_argument("--format", default="aims")
    parser.add_argument("--frac", action="store_true")
    args = parser.parse_args()

    ### Greet
    fname = args.geometry

    atoms = read(fname)

    print(f"Perfom symmetry refinement for")
    inform(atoms, symprec=args.tolerance)

    if args.prim:
        atoms = standardize_cell(atoms, to_primitve=True, symprec=args.tolerance)
        outfile = f"{args.geometry}.primitive"
    elif args.conv:
        atoms = standardize_cell(atoms, to_primitve=False, symprec=args.tolerance)
        outfile = f"{args.geometry}.conventional"
    else:
        atoms = refine_cell(atoms, symprec=args.tolerance)
        outfile = f"{args.geometry}.refined"

    write(atoms, outfile, format=args.format, spacegroup=True, scaled=args.frac)

    print(f"\nNew structure written in {args.format} format to {outfile}")


if __name__ == "__main__":
    main()
