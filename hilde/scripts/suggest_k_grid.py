from argparse import ArgumentParser as argpars
from hilde.io import read, inform
from hilde.helpers.k_grid import d2k


def main():
    """ suggest k_grid """
    parser = argpars(description="Read geometry and suggest k_grid based on density")
    parser.add_argument("geom", type=str, help="geometry input file")
    parser.add_argument("-d", "--density", type=float, default=3.5)
    parser.add_argument("--uneven", action="store_true")
    parser.add_argument("--format", default="aims")
    args = parser.parse_args()

    ### Greet
    fname = args.geom
    cell = read(fname, format=args.format)
    inform(cell, fname=fname)

    k_grid = d2k(cell, kptdensity=args.density, even=not args.uneven)

    print(f"\nSuggested k_grid for kpt-density {args.density}:")
    print(f"  k_grid ", " ".join([f"{k}" for k in k_grid]))


if __name__ == "__main__":
    main()
