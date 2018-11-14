""" Summarize output from ASE.md class (in md.log) """

from argparse import ArgumentParser
import numpy as np


def main():
    """ main routine """
    parser = ArgumentParser(description="Read md.log and make simple statistics")
    parser.add_argument("file", help="md.log input file")
    parser.add_argument("-p", "--plot", action="store_true", help="plot to pdf")
    parser.add_argument("--avg", type=int, help="running avg in plot", default=50)
    args = parser.parse_args()

    e_kin = []
    e_pot = []
    temp = []
    time = []
    with open(args.file) as f:
        for line in f:
            if line.strip() == "":
                continue
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

    if args.plot:
        plot_temperature(time, temp, args.avg)


def plot_temperature(time, temperatures, running_avg=50):
    """ Plot the nuclear temperature and save to pdf """
    import pandas as pd
    import matplotlib
    from hilde.helpers.plotting import tableau_colors as tc

    matplotlib.use('pdf')

    data = pd.Series(temperatures, time)

    ax = data.plot(x="time", y="temperatures", color=tc[0])
    ax.set_title(f"Nuclear Temperature")
    ax.legend(["Instant. Temp."])

    if running_avg > 1:
        roll = data.rolling(window=running_avg, min_periods=0).mean()
        roll.plot(x="time", y="temperatures", ax=ax, style="--", color=tc[1])
        ax.set_title(f"Nucl. Temp. with runnig avg. (window = {running_avg})")
        ax.legend(["Instant. Temp.", "Running mean"])

    ax.set_xlabel("Time [ps]")
    ax.set_ylabel("Nucl. Temperature [K]")
    fig = ax.get_figure()
    fig.savefig("temp.pdf")


if __name__ == "__main__":
    main()
