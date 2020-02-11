"""gather statistics about trajectory data"""

import pandas as pd

from vibes import keys
from vibes.green_kubo.heat_flux import (
    get_heat_flux_power_spectrum,
    get_kappa_cumulative_dataset,
)


def summary(dataset, plot=True, hann=True, **kwargs):
    """summarize heat_flux data in xarray DATASET

    Args:
        dataset(xarray.Dataset): the trajectory.dataset

    """
    assert keys.heat_flux in dataset
    assert keys.heat_flux_aux in dataset

    # total fluxes
    J1 = dataset.heat_flux.dropna(keys.time)
    J2 = dataset.heat_flux_aux.dropna(keys.time)

    # fourier transforms
    Jw1 = get_heat_flux_power_spectrum(J1, verbose=True)
    Jw2 = get_heat_flux_power_spectrum(J2, verbose=False)

    # autocorrelation
    ds = get_kappa_cumulative_dataset(dataset, hann=hann)
    ds_aux = get_kappa_cumulative_dataset(dataset, hann=hann, aux=True, verbose=False)

    # time resolved
    d = {
        "Jcorr": ds[keys.hfacf_scalar],
        "kappa": ds[keys.kappa_cumulative_scalar],
        "kappa_aux": ds_aux[keys.kappa_cumulative_scalar],
    }
    df_time = pd.DataFrame(d, index=ds.time)

    # freq resolved
    d = {"Jspec": Jw1.sum(axis=1), "Jspec_aux": Jw2.sum(axis=1)}
    df_freq = pd.DataFrame(d, index=Jw1.omega)

    return (df_time, df_freq)


def plot_summary(df_time, df_freq, avg=50, logx=True, xlim=None):
    """plot a summary of the data in DATAFRAME"""
    import matplotlib
    from matplotlib import pyplot as plt
    from vibes.helpers.plotting import tableau_colors as tc

    matplotlib.use("pdf")

    try:
        import seaborn as sns

        sns.set_style("whitegrid")
    except ModuleNotFoundError:
        pass

    # settings for the immediate plot
    alpha = 0.5
    color = tc[3]
    plot_kw = {
        "alpha": alpha,
        "linewidth": 0.5,
        "label": "",
        "color": color,
        "marker": ".",
    }
    avg_kw = {"linewidth": 3, "color": "k"}
    fig_kw = {
        "figsize": (11.69, 8.27),
        "gridspec_kw": {"width_ratios": [3, 2]},
        "sharex": "col",
    }
    kw_roll = {"window": avg, "min_periods": 0, "center": True}

    df_time.index /= 1000

    jc = df_time.Jcorr / df_time.Jcorr.iloc[0]
    kc = df_time.kappa
    kc_aux = df_time.kappa_aux
    js = df_freq.Jspec
    js_aux = df_freq.Jspec_aux

    fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(nrows=2, ncols=2, **fig_kw)

    # plot
    # HFACF
    jc.plot(ax=ax1, **plot_kw)
    jc.rolling(**kw_roll).mean().plot(ax=ax1, **avg_kw)

    # Kappa
    kc_aux.plot(ax=ax2, **plot_kw)
    kc_aux.rolling(**kw_roll).mean().plot(ax=ax2, **avg_kw)
    kc.plot(ax=ax2, **{**avg_kw, "color": tc[1], "linewidth": 2})

    # plot spectra
    js.plot(ax=ax3, **plot_kw)
    js_aux.plot(ax=ax4, **plot_kw)
    js.rolling(**kw_roll).mean().plot(ax=ax3, **avg_kw)
    js_aux.rolling(**kw_roll).mean().plot(ax=ax4, **avg_kw)

    if xlim:
        ax2.set_xlim((-1, xlim))
    else:
        ax2.set_xlim((-1, kc.index.max()))
    ax2.set_ylim((0, 1.1 * kc.max()))
    if logx:
        ax2.set_xscale("log")
        ax2.set_xlim(kc.index[2], kc.index.max())

    linthreshy = 0.1 * js_aux.max()
    ax4.set_yscale("symlog", linthreshy=linthreshy)
    ax4.axhline(linthreshy, ls="--", color="purple")
    ax4.set_xlim((-1, js_aux.index.max()))

    # titles and labels
    fig.suptitle(f"Heat Flux Overview (window: {avg})")
    ax1.set_title(r"$\langle J (t) J(0) \rangle ~/~ \langle J (0) J(0) \rangle$")
    ax2.set_title(r"$\kappa (t)$ (W/mK)")
    ax2.set_xlabel("Time $t$ (ps)")
    ax3.set_title(r"$\vert J (\omega)\vert^2$")
    ax4.set_xlabel("Omega (THz)")
    ax4.set_title(r"$\vert J_\mathrm{aux} (\omega)\vert^2$")

    return fig
