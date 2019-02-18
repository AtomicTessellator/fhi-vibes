""" a wrapper for TDEP """
from ase.io.aims import read_aims
from subprocess import run
from pathlib import Path

from hilde.helpers.paths import cwd
from hilde.helpers import Timer
from hilde.phonopy.postprocess import extract_results
from hilde.trajectory import reader

def convert_phonopy_to_dep(
    ph, workdir="tdep", logfile="convert_phonopy_to_dep.log"
):
    command = ["convert_phonopy_to_forceconstant", "--truncate"]
    with cwd(workdir, mkdir=True), open(logfile, "w") as file:
        extract_results(ph, tdep=True)
        run(command, stdout=file)
        outfile = Path("outfile.converted_forceconstant")
        infile = Path("infile.forceconstant")
        if infile.exists():
            proceed = input(f"Symlink {infile} exists. Proceed? (y/n) ")
            if proceed.lower() == "y":
                infile.unlink()
            else:
                print(".. Symlink NOT created.")
                return

        infile.symlink_to(outfile)
        print(f".. Symlink {infile} created.")

def generate_cannonical_configurations(
     ph=None, workdir="tdep", temperature=300, n_sample=5, quantum=False, logfile="canon_conf.log"
):
    if ph:
        convert_phonopy_to_dep(ph, workdir)
    command = ["canonical_configuration"]
    if quantum:
        command.append(f"--quantum")
    command.extend("-of 4".split())
    command.extend(f"-n {n_sample}".split())
    command.extend(f"-t {temperature}".split())
    with cwd(workdir, mkdir=True), open(logfile, "w") as file:
        run(command, stdout=file)
    outfiles = Path(workdir).glob("aims_conf*")
    return [read_aims(of) for of in outfiles]

def extract_forceconstants_from_trajectory(
    trajectory_file, workdir="tdep", rc2=10, remapped=True, logfile="fc.log", **kwargs
):
    trajectory = reader(trajectory_file)
    if "skip" in kwargs:
        skip = kwargs["skip"]
    else:
        skip = 0
    trajectory.to_tdep(folder=workdir, skip=0)
    extract_forceconstants(workdir, rc2, remapped, logfile, **kwargs)

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
