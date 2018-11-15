#!/usr/bin/env python2.7

# USAGE:  ./get_relaxation_info.py  aims.out (aims.out.2 ...)
#
# Revision 2018/08: FK

from argparse import ArgumentParser as argpars

parser = argpars(description="Summarize the relaxation path")
parser.add_argument("aimsouts", type=str, nargs="+", help="aims output files")
args = parser.parse_args()

# Find the optimizer type
def get_optimizer(f):
    line = next(l for l in f if "Geometry relaxation:" in l)
    if "Textbook BFGS" in line:
        return 1
    elif "TRM" in line:
        return 2


# find energy
def get_energy(f):
    try:
        line = next(l for l in f if "Total energy uncorrected" in l)
        total_energy = float(line.split()[5])
    except StopIteration:
        exit()
    line = next(l for l in f if "Electronic free energy" in l)
    free_energy = float(line.split()[5])
    return total_energy, free_energy


# get max_force
def get_forces(f):
    line = next(l for l in f if "Maximum force component" in l)
    return float(line.split()[4])


# parse info of one step
def parser(f, n_init=0, optimizer=2):
    n_rel = n_init
    converged = 0
    abort = 0
    while not converged and not abort:
        n_rel += 1
        status = 0
        energy, free_energy = get_energy(f)
        max_force = get_forces(f)
        for line in f:
            if "Present geometry is converged." in line:
                converged = 1
                break
            elif "Advancing" in line:
                pass
            elif "Aborting optimization" in line:
                abort = 1
            elif "Counterproductive step -> revert!" in line:
                status = 1
            elif "Optimizer is stuck" in line:
                status = 2
            #            elif '**' in line:
            #                status = 3
            elif "Finished advancing geometry" in line:
                break
            elif "Updated atomic structure" in line:
                break
        yield n_rel, energy, free_energy, max_force, status, converged, abort


def print_status(n_rel, energy, de, free_energy, df, max_force, status_string):
    print(
        "{:5d}   {:16.8f} {:14.6f}   {:16.8f} {:14.6f} {:10.6f}  {}".format(
            n_rel, energy, de, free_energy, df, max_force, status_string
        )
    )


def main():
    init, n_rel, converged, abort = 4 * (None,)
    status_string = [
        "",
        "rejected.",
        "rejected: force <-> energy inconsistency?",
        "stuck.",
    ]

    # Run
    print(
        "\n# Step Total energy [eV]   E-E(1) [meV]   Free energy [eV]   F-F(1)"
        + " [meV]   max. force [eV/AA]\n"
    )

    for infile in args.aimsouts:
        with open(infile) as f:
            # Check optimizer
            optimizer = get_optimizer(f)
            ps = parser(f, n_init=n_rel or 0, optimizer=optimizer)
            for (n_rel, energy, free_energy, max_force, status, converged, abort) in ps:
                if not init:
                    first_energy, first_free_energy = energy, free_energy
                    init = 1
                print_status(
                    n_rel,
                    energy,
                    1000 * (energy - first_energy),
                    free_energy,
                    1000 * (free_energy - first_free_energy),
                    max_force,
                    status_string[status],
                )

    if converged:
        print("--> converged.")
    if abort:
        print("*--> aborted, too many steps.")


if __name__ == "__main__":
    main()
