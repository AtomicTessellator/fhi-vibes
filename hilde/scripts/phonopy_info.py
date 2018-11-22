""" Summarize output from ASE.md class (in md.log) """

from argparse import ArgumentParser
from pathlib import Path


def preprocess(args):
    import numpy as np
    from ase.io import read
    from hilde.settings import Settings
    import hilde.phonopy.wrapper as ph

    atoms = read(args.infile, format=args.format)

    _, _, scs_ref = ph.preprocess(atoms, supercell_matrix=1)

    if args.dim is not None:
        phonon, sc, scs = ph.preprocess(atoms, supercell_matrix=args.dim)
    else:
        settings = Settings(args.config_file)
        phonon, sc, scs = ph.preprocess(atoms, **settings.phonopy)
        print("Phonopy settings:")
        settings.print()

    sc_str = np.array2string(phonon.get_supercell_matrix().flatten(), separator=", ")
    print("Phonopy Information")
    print(f"  Supercell matrix:        {sc_str}")
    print(f"  Number of atoms in SC:   {len(sc)}")
    print(f"  Number of displacements: {len(scs)} ({len(scs_ref)})")


def postprocess(args):
    from hilde.helpers.pickletools import pread
    from hilde.phonopy.wrapper import plot_bandstructure

    phonon = pread(args.infile)

    plot_bandstructure(phonon, file="bandstructure.pdf")


def main():
    """ main routine """
    parser = ArgumentParser(description="information about phonopy task")
    parser.add_argument("infile", help="primitive structure or pickled phonopy")
    parser.add_argument("--dim", type=int, nargs="*", default=None)
    parser.add_argument("--config_file", default="phonopy.cfg")
    parser.add_argument("--format", default="aims")
    args = parser.parse_args()

    suffix = Path(args.infile).suffix
    if suffix == ".py":
        preprocess(args)

    elif suffix == ".pick" or suffix == ".gz":
        postprocess(args)


if __name__ == "__main__":
    main()
