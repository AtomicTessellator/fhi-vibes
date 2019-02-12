""" a wrapper for TDEP """

from subprocess import run
from pathlib import Path

from hilde.helpers.paths import cwd
from hilde.helpers import Timer


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
                print(f'.. use gnuplot to plot dispersion to pdf')
                command = "gnuplot -p outfile.dispersion_relations.gnuplot_pdf".split()
                run(command, stdout=file)

    timer()
