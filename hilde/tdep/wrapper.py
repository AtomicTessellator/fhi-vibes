""" a wrapper for TDEP """

from subprocess import run
from pathlib import Path

import numpy as np

from hilde.helpers.paths import cwd
from hilde.helpers import Timer


def parse_tdep_forceconstant(fname="infile.forceconstant_remapped", remapped=True):
    """ parse the remapped forceconstants from TDEP """
    timer = Timer()

    if not remapped:
        raise RuntimeError("Only remapped forceconstants can be parsed currently.")

    print(f"Parse force constants from\n  {fname}")
    with open(fname) as fo:
        n_atoms = int(next(fo).split()[0])
        cutoff = float(next(fo).split()[0])

        print(f".. Number of atoms:   {n_atoms}")
        print(rf".. Real space cutoff: {cutoff:.3f} \AA")

        force_constants = np.zeros([n_atoms, 3, n_atoms, 3])

        for i1 in range(n_atoms):
            n_neighbors = int(next(fo).split()[0])
            for _ in range(n_neighbors):
                i2 = int(next(fo).split()[0]) - 1
                # skip the lattice point
                _ = np.array(next(fo).split(), dtype=float)
                phi = np.array([next(fo).split() for _ in range(3)], dtype=float)

                force_constants[i1, :, i2, :] = phi

    timer()

    return force_constants.reshape(2 * (3 * n_atoms,))


def extract_forceconstants(
    workdir="tdep", rc2=10, remapped=True, logfile="fc.log", **kwargs
):
    """ run tdep's extract_forceconstants in the working directory """

    timer = Timer()

    print(f"Extract force constants with TDEP from input files in\n  {workdir}")

    command = ["extract_forceconstants", "--verbose"]

    command.extend(f"-rc2 {rc2}".split())

    if remapped:
        command.append("--printfc2remapped")

    with cwd(workdir), open(logfile, "w") as file:
        run(command, stdout=file)
        timer()

        # create the symlink of force constants
        print(f".. Create symlink to forceconstant file")
        outfile = Path("outfile.forceconstant")
        infile = Path("infile" + outfile.suffix)

        if infile.exists():
            proceed = input(f"Symlink {infile} exists. Proceed? (y/n) ")
            if proceed.lower() == "y":
                infile.unlink()
            else:
                print(".. Symlink NOT created.")
                return

        infile.symlink_to(outfile)
        print(f".. Symlink {infile} created.")


def phonon_dispersion_relations(workdir="tdep", gnuplot=True, logfile="dispersion.log"):
    """ run tdep's phonon_dispersion_relations in working directory """

    timer = Timer(f"Run TDEP phonon_dispersion_relations in {workdir}")

    with cwd(workdir):
        # check if input files are present
        for file in ("forceconstant", "ucposcar", "ssposcar"):
            path = Path("infile." + file)
            if not path.exists():
                raise FileNotFoundError(f"{path} missing in ./{workdir}.")

        # plot if ipnut files are present
        command = "phonon_dispersion_relations -p".split()

        with open(logfile, "w") as file:
            run(command, stdout=file)

            if gnuplot:
                print(f".. use gnuplot to plot dispersion to pdf")
                command = "gnuplot -p outfile.dispersion_relations.gnuplot_pdf".split()
                run(command, stdout=file)

    timer()
