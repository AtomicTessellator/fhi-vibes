import os
from argparse import ArgumentParser as argpars
from hilde.io import read, write
from hilde.helpers.structure import clean_atoms
from hilde.structure.io import inform


def main():
    """ print geometry information """
    parser = argpars(description="Read geometry and print some symmetry info")
    parser.add_argument("geom", type=str, help="geometry input file")
    parser.add_argument("--align", action="store_true", help="align the lattice")
    parser.add_argument("--format", default="aims")
    args = parser.parse_args()

    atoms = read(args.geom, format=args.format)
    inform(atoms)

    atoms = clean_atoms(atoms, align=args.align)

    outfile = f"{args.geom}.cleaned"
    write(atoms, outfile, format=args.format)

    print(f"\nCleaned geometry written to {outfile}")


if __name__ == "__main__":
    main()
