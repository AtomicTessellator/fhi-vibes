""" Summarize output from ASE.md class (in md.log) """

from pathlib import Path
from argparse import ArgumentParser
import numpy as np
from hilde.trajectory import reader


def parse_log(filename):
    """ parse the ASE logfile, typically md.log """
    e_kin = []
    e_pot = []
    temp = []
    time = []
    with open(filename) as f:
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

    return e_kin, e_pot, temp, time


def main():
    """ main routine """
    parser = ArgumentParser(description="Read md.log and make simple statistics")
    parser.add_argument("file", help="md.log or trajectory.yaml input file")
    parser.add_argument("-p", "--plot", action="store_true", help="plot to pdf")
    parser.add_argument("--avg", type=int, help="running avg in plot", default=50)
    parser.add_argument("-v", "--verbose", action="store_true", help="give more info")
    args = parser.parse_args()

    infile = Path(args.file)

    if "yaml" in infile.suffix:
        trajectory = reader(infile)[1:]
        e_kin = [atoms.get_kinetic_energy() for atoms in trajectory]
        e_pot = [atoms.get_potential_energy() for atoms in trajectory]
        temp = [atoms.get_temperature() for atoms in trajectory]
        # sum up time steps, necessary if dt changes during md
        time = np.cumsum([atoms.info["dt"] * 0.001 for atoms in trajectory])
    elif "log" in infile.suffix:
        e_kin, e_pot, temp, time = parse_log(infile)

    if args.verbose:
        n_range = min(20, len(trajectory) - 2)
        print(f"step   time   temperature")
        for ii in range(n_range):
            print(
                f'{trajectory[-n_range+ii].info["nsteps"]:5d} '
                f"{time[-n_range+ii]:10.4f} "
                f"{temp[-n_range+ii]:10.4f} "
            )

    # divide into three thirds
    len_3 = len(temp) // 3
    temp_1 = temp[:len_3]
    temp_2 = temp[len_3 : 2 * len_3]
    temp_3 = temp[2 * len_3 :]

    print(f"Simulation time:         {time[-1] - time[0]:.4f} ps ({len(time)} steps)")
    print(f"Temperature:             {np.mean(temp):.2f} +/- {np.std(temp):.2f}K")
    print(f"Temperature (1st 1/3):   {np.mean(temp_1):.2f} +/- {np.std(temp_1):.2f}K")
    print(f"Temperature (2nd 1/3):   {np.mean(temp_2):.2f} +/- {np.std(temp_2):.2f}K")
    print(f"Temperature (3rd 1/3):   {np.mean(temp_3):.2f} +/- {np.std(temp_3):.2f}K")
    print(f"Kinetic energy:          {np.mean(e_kin):.2f} +/- {np.std(e_kin):.2f}eV")
    print(f"Potential energy:        {np.mean(e_pot):.2f} +/- {np.std(e_pot):.2f}eV")

    if args.plot:
        plot_temperature(time, temp, args.avg)


def plot_temperature(time, temperatures, running_avg=50):
    """ Plot the nuclear temperature and save to pdf """
    import pandas as pd
    import matplotlib
    from hilde.helpers.plotting import tableau_colors as tc

    matplotlib.use("pdf")

    data = pd.Series(temperatures, time)

    temp_2 = data.iloc[len(data) // 2 :]

    ax = data.plot(
        x="time",
        y="temperatures",
        color=tc[0],
        title=f"Nuclear Temperature",
        legend=["Instant. Temp."],
    )
    # ax.legend(["Instant. Temp."])

    if running_avg > 1:
        roll = data.rolling(window=running_avg, min_periods=0).mean()
        roll.plot(
            x="time",
            y="temperatures",
            ax=ax,
            style="--",
            color=tc[1],
            title=(
                f"Nucl. Temp. with runnig avg. (window = {running_avg})"
                f"\nTemperature (last 1/2): {temp_2.mean():.2f} +/- {temp_2.std():.2f}K"
            ),
            legend=["Instant. Temp.", "Running mean"],
        )

    ax.set_xlabel("Time [ps]")
    ax.set_ylabel("Nucl. Temperature [K]")
    fig = ax.get_figure()
    fig.savefig("temp.pdf")


if __name__ == "__main__":
    main()
