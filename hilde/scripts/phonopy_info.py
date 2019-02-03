""" Summarize output from ASE.md class (in md.log) """

import numpy as np
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
    print(f"  Superlattice:")
    for latvec in sc.cell:
        lv_str = "{:-6.2f} {:-6.2f} {:-6.2f}".format(*latvec)
        print(f"                         {lv_str}")
    print(f"  Number of atoms in SC:   {len(sc)}")
    print(f"  Number of displacements: {len(scs)} ({len(scs_ref)})")


def postprocess(args):
    from hilde.helpers.pickle import pread
    from hilde.phonopy.wrapper import plot_bandstructure, plot_bandstructure_and_dos
    from hilde.phonopy.wrapper import summarize_bandstructure, get_force_constants

    phonon = pread(args.infile)

    plot_bandstructure(phonon, file="bandstructure.pdf")
    if args.dos == "t":
        plot_bandstructure_and_dos(phonon, file="bands_and_dos.pdf")
    elif args.dos == "p":
        plot_bandstructure_and_dos(phonon, partial=True, file="bands_and_pdos.pdf")
    summarize_bandstructure(phonon, fp_file=args.fp_file)

    force_costants = get_force_constants(phonon)
    filname = "force_constants.dat"
    np.savetxt(filname, force_costants)
    print(f"Force constants saved to {filname}.")


def main():
    """ main routine """
    parser = ArgumentParser(description="information about phonopy task")
    parser.add_argument("infile", help="primitive structure or pickled phonopy")
    parser.add_argument("--dim", type=int, nargs="*", default=None)
    parser.add_argument("--config_file", default="settings.in")
    parser.add_argument("--format", default="aims")
    parser.add_argument("--fp_file", default=None, help="File to store the fingerprint")
    parser.add_argument(
        "--dos", nargs="?", const="t", default=None, help="plot dos as well"
    )
    args = parser.parse_args()

    suffix = Path(args.infile).suffix
    if suffix == ".in":
        preprocess(args)

    elif suffix == ".pick" or suffix == ".gz":
        postprocess(args)
    else:
        print("*** Nothing happened.")


if __name__ == "__main__":
    main()
