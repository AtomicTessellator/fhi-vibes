"""gather statistics about trajectory data"""
import numpy as np

from vibes import keys

# from vibes.correlation import get_autocorrelation_exponential
from vibes.helpers.xarray import xtrace


def plot_summary(dataset, avg=50, logx=True, xlim=None):
    """plot a summary of the data in DATASET"""
    import matplotlib
    from matplotlib import pyplot as plt

    # from vibes.helpers.plotting import tableau_colors as tc

    matplotlib.use("pdf")

    try:
        import seaborn as sns

        sns.set_style("white")
        sns.set_palette("colorblind")
    except ModuleNotFoundError:
        pass

    dataset = dataset.mean(dim=keys.trajectory)

    # settings for the immediate plot
    alpha = 0.5
    color = "k"  # tc[3]
    kw1 = {"color": color, "alpha": alpha}
    kw2 = {"color": color, "linewidth": 2}
    kw_roll = {"min_periods": 5, "center": True}
    gridspec = {"width_ratios": [3, 2]}
    fig_kw = {"figsize": (11.69, 8.27), "gridspec_kw": gridspec, "sharex": "col"}

    fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(nrows=2, ncols=2, **fig_kw)

    dataset[keys.time] = dataset[keys.time] / 1000
    time = dataset[keys.time]

    # plot 3 hf autocorrelation functions
    j_corr = dataset[keys.heat_flux_autocorrelation]
    for ii in range(3):
        jc = j_corr[:, ii, ii]
        jc.rolling({keys.time: 5}, **kw_roll).mean().to_series().plot(ax=ax11, lw=0.5)
    jc = xtrace(j_corr) / 3
    jc.rolling({keys.time: 5}, **kw_roll).mean().to_series().plot(ax=ax11, **kw2)
    ax11.set_ylim([j_corr.min(), j_corr[min(2, len(j_corr)) - 1].max() * 1.1])

    # plot 3 kappas
    kappa = dataset[keys.kappa_cumulative]
    # t_avalanche = dataset[keys.time_avalanche] / 1000
    for ii in range(3):
        kc = kappa[:, ii, ii]  # .sel({keys.time: slice(0, t)})
        kc.rolling({keys.time: 5}, **kw_roll).mean().to_series().plot(ax=ax21, **kw1)
        ax21.axhline(dataset.kappa[ii, ii], c="k", alpha=alpha, lw=0.2)
        # ax21.axvline(t_avalanche[ii, ii], c="k", alpha=alpha)
    k_diag = np.diag(dataset.kappa)
    mean = k_diag.mean()
    dev = k_diag.std() / 3 ** 0.5
    ax21.axhline(mean, c="k", lw=2)
    ax21.fill_between(time, mean - dev, mean + dev, color="k", alpha=0.3)
    kc = xtrace(kappa) / 3
    kc.rolling({keys.time: 5}, **kw_roll).mean().to_series().plot(ax=ax21, **kw2)

    ax21.set_xlim([-1, xlim or time.max()])
    # ax21.set_xscale("log")

    # plot power spectrum
    pow = dataset[keys.heat_flux_power_spectrum]
    for ii in range(3):
        jw = pow[:, ii, ii]
        jw.to_series().plot(ax=ax12, **kw1)
        adjust_y(ax12, jw)

    if keys.heat_flux_total_power_spectrum in dataset:
        pow = dataset[keys.heat_flux_total_power_spectrum]
        for ii in range(3):
            jw = pow[:, ii, ii].real
            jw.to_series().plot(ax=ax22, **kw1)
            adjust_y(ax22, jw)
        ax22.set_yscale("log")
        ax22.set_title(r"$\vert J_\mathrm{aux} (\omega)\vert^2$")
    else:
        ax22.set_title(r"$J_\mathrm{aux} (\omega)$ missing.")

    # titles and labels
    fig.suptitle(f"Heat Flux Overview (window: {avg})")
    ax11.set_title(r"$\langle J (t) J(0) \rangle ~/~ \langle J (0) J(0) \rangle$")
    ax21.set_title(r"$\kappa (t)$ (W/mK)")
    ax21.set_xlabel("Time $t$ (ps)")
    ax12.set_title(r"$\vert J (\omega)\vert^2$")
    ax22.set_xlabel("Omega (THz)")
    return fig


def adjust_y(ax, data, positive=True):
    y0, y1 = ax.get_ylim()
    y0 = abs(min(y0, np.min(abs(data)) * 0.9))
    y1 = max(y1, 1.2 * np.max(abs(data)))
    ax.set_ylim([y0, y1])
