""" Summarize output from ASE.md class (in md.log) """

from argparse import ArgumentParser
import numpy as np


def main():
    parser = ArgumentParser(description="Read md.log and make simple statistics")
    parser.add_argument("file", help="md.log input file")
    args = parser.parse_args()

    e_kin = []
    e_pot = []
    temp = []
    time = []
    with open(args.file) as f:
        for line in f:
            if "Time" in line:
                continue
            t, _, ep, ek, T = (float(l) for l in line.split())
            time.append(t)
            temp.append(T)
            e_kin.append(ek)
            e_pot.append(ep)

    print(f"Simulation time:        {time[-1] - time[0]:.4f}ps")
    print(f"Temperature:            {np.mean(temp):.2f} +/- {np.std(temp):.2f}K")
    print(f"Kinetic energy:         {np.mean(e_kin):.2f} +/- {np.std(e_kin):.2f}eV")
    print(f"Potential energy:       {np.mean(e_pot):.2f} +/- {np.std(e_pot):.2f}eV")


if __name__ == "__main__":
    main()
