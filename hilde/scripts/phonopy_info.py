""" Summarize output from ASE.md class (in md.log) """

import json
from argparse import ArgumentParser
import numpy as np
from ase.io import read
from hilde.settings import Settings
import hilde.phonopy.phono as ph


def main():
    """ main routine """
    parser = ArgumentParser(description="convert yaml file to json")
    parser.add_argument("geometry", help="primitive structure")
    parser.add_argument("--dim", type=int, nargs="*", default=None)
    parser.add_argument("--config_file", default="phonopy.cfg")
    parser.add_argument("--format", default="aims")
    args = parser.parse_args()

    atoms = read(args.geometry, format=args.format)

    _, _, scs_ref = ph.preprocess(atoms, supercell_matrix=1)

    if args.dim is not None:
        phonon, _, scs = ph.preprocess(atoms, supercell_matrix=args.dim)
    else:
        settings = Settings(args.config_file)
        phonon, _, scs = ph.preprocess(atoms, **settings.phonopy)
        print("Phonopy settings:")
        settings.print()

    sc_str = np.array2string(phonon.get_supercell_matrix().flatten(), separator=", ")
    print("Phonopy Information")
    print(f"  Supercell matrix:        {sc_str}")
    print(f"  Number of displacements: {len(scs)} ({len(scs_ref)})")


if __name__ == "__main__":
    main()
