"""plot trajectory data"""
from ase import units

from vibes import keys


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

    from vibes.helpers.plotting import tableau_colors as tc

    # plot temperatures
    temp = dataframe[keys.temperature]
    e_kin = dataframe[keys.energy_kinetic]
    e_pot = dataframe[keys.energy_potential]

    e_pot -= e_pot.min()

    # settings for the immediate plot
    plot_kw = {"alpha": 0.4, "linewidth": 1.0, "label": ""}
    avg_kw = {"linewidth": 1.5}
    fig_kw = {
        "figsize": (8.27, 11.69),
        "gridspec_kw": {"height_ratios": [1, 1, 0.5, 0.5]},
        "sharex": True,
    }
    kw_roll = {"window": avg, "min_periods": 0, "center": True}

    # pressure: make sure there is enough data, otherwise don't bother plotting
    p = dataframe.pressure / units.GPa
    p_int = p.interpolate("akima")
    p = p.dropna()
    if len(p) > 3:
        fig, (ax, ax2, ax3, ax4) = plt.subplots(nrows=4, **fig_kw)
        # fig = plt.figure(**fig_kw)
        # gs = gridspec.GridSpec(nrows=6, ncols=1)
        # ax = fig.add_subplot(gs[0:2])
        # ax2 = fig.add_subplot(gs[2:4])
        # ax3 = fig.add_subplot(gs[4])
        # ax4 = fig.add_subplot(gs[5])
        # fig, (ax, ax2, ax3) = plt.subplots(nrows=3, **fig_kw)
        p.plot(
            ax=ax3, **{**plot_kw, "label": "p", "linewidth": 0}, color=tc[3], marker="x"
        )
        p_int.plot(ax=ax3, **{**plot_kw, "label": "p (akima)"}, color=tc[3])
        p_int.expanding().mean().plot(ax=ax3, label="avg. p", **avg_kw, color=tc[3])
        ax3.axhline(0, linewidth=0.75)
        ax3.legend(loc=4)
        ax3.set_title("Pressure")
        ax3.set_ylabel("Pressure [GPa]")

    else:
        fig_kw["gridspec_kw"] = {"height_ratios": [1, 1, 1]}
        fig, (ax, ax2, ax4) = plt.subplots(nrows=3, **fig_kw)

    temp.plot(color=tc[3], title="Nuclear Temperature", ax=ax, **plot_kw)

    if natoms:
        e_temp = (e_kin + e_pot) / natoms / 3 / units.kB
        e_temp.plot(color=tc[1], ax=ax, **plot_kw)
        e_temp.rolling(**kw_roll).mean().plot(
            color=tc[1], label=f"E_tot", ax=ax, **avg_kw
        )

    roll = temp.rolling(**kw_roll).mean()
    roll.plot(color=tc[3], label=f"T_nucl", ax=ax, **avg_kw)
    ax.axhline(temp.mean(), linewidth=0.75)

    # exp = data.iloc[min(len(data) // 2, args.avg) :].expanding().mean()
    # exp.plot( color=tc[5], label=f"Expanding mean ({args.avg}", ax=ax)

    ax.set_xlabel("Time [ps]")
    ax.set_ylabel("Nucl. Temperature [K]")
    ax.legend(loc=4)

    # fig.savefig("temp.pdf")

    # plot energies in one plot
    # fig, ax = plt.subplots()

    e_tot = e_pot + e_kin
    # e_dif = e_pot #e_kin - e_pot

    e_tot.plot(color=tc[0], title="Energy", ax=ax2, **plot_kw)
    roll = e_tot.rolling(**kw_roll).mean()
    roll.plot(color=tc[0], ax=ax2, label="E_tot", **avg_kw)

    e_kin.plot(color=tc[3], ax=ax2, **plot_kw)
    roll = e_kin.rolling(**kw_roll).mean()
    roll.plot(color=tc[3], ax=ax2, label="E_kin", **avg_kw)

    ax2.axhline(0, linewidth=0.75)
    e_pot.plot(color=tc[1], ax=ax2, **plot_kw)
    exp = e_pot.rolling(**kw_roll).mean()
    exp.plot(color=tc[1], ax=ax2, label="E_pot", **avg_kw)

    ax2.legend(loc=4)
    ax2.set_ylabel("Energy [eV]")

    # displacements
    dr = dataframe.dr_mean
    dr_std = dataframe.dr_std
    dr.plot(color=tc[4], ax=ax4, **plot_kw)
    roll = dr.rolling(**kw_roll).mean()
    roll.plot(color=tc[4], ax=ax4, label="mean(|dR|)", **avg_kw)

    dr_std.plot(color=tc[5], ax=ax4, **plot_kw)
    roll = dr_std.rolling(**kw_roll).mean()
    roll.plot(color=tc[5], ax=ax4, label="std(|dR|)", **avg_kw)
    # ax4.fill_between(dr.index, dr + dr_std, dr - dr_std, color=tc[4], alpha=0.3)
    ax4.set_ylabel("Displacement [Å]")
    ax4.set_title("Displacement")
    ax4.axhline(0)
    ax4.legend(loc=4)

    ax4.set_xlabel("Time [ps]")

    # fig.tight_layout()
    fname = "md_summary.pdf"
    fig.savefig(fname, bbox_inches="tight")
    print(f".. summary plotted to {fname}")
