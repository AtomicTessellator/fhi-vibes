from hilde.parsers import read_structure
from argparse import ArgumentParser as argpars
import os


def main():
    parser = argpars(description="Read geometry and make primitive + conventional unit")
    parser.add_argument("geom", type=str, help="geometry input file")
    parser.add_argument(
        "-t",
        "--tolerance",
        type=float,
        default=1e-5,
        help="symmetry tolerance (symprec)",
    )
    parser.add_argument(
        "--aimsout",
        nargs="?",
        type=str,
        default=None,
        const="-1",
        help="use aims output file",
    )
    parser.add_argument("--dry", action="store_true", help="only show symmetry")
    parser.add_argument(
        "-c", "--cartesian", action="store_true", help="write cartesian positions"
    )
    parser.add_argument("--conv", action="store_true", help="store conventional cell")
    args = parser.parse_args()

    ### Greet
    fname = args.geom
    cell = read_structure(fname, symprec=args.tolerance)
    print(f"Perfom symmetry refinement on")
    print(f"  input geometry:  {cell.sysname}")
    print(f"  from:            {fname}")

    # Make primitive standardized:
    prim_cell = cell.get_primitive_standardized()
    # Make conventional standardized:
    conv_cell = cell.get_conventional_standardized()

    if args.dry:
        exit()

    scaled_pos = not args.cartesian
    if args.conv:
        conv_cell.write(fname + ".conv", scaled=scaled_pos)
        print(f"  conventional cell written to: {fname}.conv")
    else:
        prim_cell.write(fname + ".prim", scaled=scaled_pos)
        print(f"  primitive cell written to:    {fname}.prim")
