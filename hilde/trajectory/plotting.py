"""plot trajectory data"""
from ase import units


def plot_summary(dataframe, avg=50, natoms=None):

    """plot a summary of the data in DATAFRAME

    Args:
        dataframe (pandas.Dataframe): MD data
        avg (int): window size for averaging
        natoms (int): number of atoms

    Returns:
        an A4 plot in `md_summary.pdf`
    """
    import matplotlib

    matplotlib.use("pdf")

    from matplotlib import pyplot as plt

    from hilde.helpers.plotting import tableau_colors as tc

    # plot temperatures
    temp = dataframe.temperature
    e_kin = dataframe.kinetic_energy
    e_pot = dataframe.potential_energy
    e_pot -= e_pot.min()

    # settings for the immediate plot
    plot_kw = {"alpha": 0.4, "linewidth": 1.0, "label": ""}
    avg_kw = {"linewidth": 1.5}
    fig_kw = {"sharex": True, "figsize": (8.27, 11.69)}

    # pressure: make sure there is enough data, otherwise don't bother plotting
    p = dataframe.pressure / units.GPa
    p_int = p.interpolate("akima")
    p = p.dropna()
    if len(p) > 3:
        fig, (ax, ax2, ax3) = plt.subplots(nrows=3, **fig_kw)
        p.plot(
            ax=ax3, **{**plot_kw, "label": "p", "linewidth": 0}, color=tc[3], marker="x"
        )
        p_int.plot(ax=ax3, **{**plot_kw, "label": "p (akima)"}, color=tc[3])
        p_int.expanding().mean().plot(ax=ax3, label="avg. p", **avg_kw, color=tc[3])
        ax3.axhline(0, linewidth=0.75)
        ax3.legend()
        ax3.set_title("Pressure")
        ax3.set_xlabel("Time [ps]")
        ax3.set_ylabel("Pressure [GPa]")
    else:
        fig, (ax, ax2) = plt.subplots(nrows=2, **fig_kw)

    temp.plot(color=tc[3], title="Nuclear Temperature", ax=ax, **plot_kw)

    if natoms:
        e_temp = (e_kin + e_pot) / natoms / 3 / units.kB
        e_temp.plot(color=tc[1], ax=ax, **plot_kw)
        e_temp.rolling(window=avg, min_periods=0).mean().plot(
            color=tc[1], label=f"E_tot", ax=ax, **avg_kw
        )

    roll = temp.rolling(window=avg, min_periods=0).mean()
    roll.plot(color=tc[3], label=f"T_nucl", ax=ax, **avg_kw)
    ax.axhline(temp.mean(), linewidth=0.75)

    # exp = data.iloc[min(len(data) // 2, args.avg) :].expanding().mean()
    # exp.plot( color=tc[5], label=f"Expanding mean ({args.avg}", ax=ax)

    ax.set_xlabel("Time [ps]")
    ax.set_ylabel("Nucl. Temperature [K]")
    ax.legend()

    # fig.savefig("temp.pdf")

    # plot energies in one plot
    # fig, ax = plt.subplots()

    e_tot = e_pot + e_kin
    e_dif = e_pot - e_kin

    e_tot.plot(color=tc[0], title="Energy", ax=ax2, **plot_kw)
    roll = e_tot.rolling(window=avg, min_periods=0).mean()
    roll.plot(color=tc[0], ax=ax2, label="E_tot", **avg_kw)

    e_pot.plot(color=tc[3], ax=ax2, **plot_kw)
    roll = e_pot.rolling(window=avg, min_periods=0).mean()
    roll.plot(color=tc[3], ax=ax2, label="E_pot", **avg_kw)

    ax2.axhline(0, linewidth=0.75)
    e_dif.plot(color=tc[1], ax=ax2, **plot_kw)
    exp = e_dif.rolling(min_periods=0, window=avg).mean()
    exp.plot(color=tc[1], ax=ax2, label="E_pot - E_kin", **avg_kw)

    ax2.legend()
    ax2.set_ylabel("Energy [eV]")

    # fig.tight_layout()
    fname = "md_summary.pdf"
    fig.savefig(fname, bbox_inches="tight")
    print(f".. summary plotted to {fname}")
