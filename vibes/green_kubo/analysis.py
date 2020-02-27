"""gather statistics about trajectory data"""
import pandas as pd

from vibes import keys
from vibes.correlation import get_autocorrelation_exponential
from vibes.fourier import get_fourier_transformed
from vibes.helpers import xtrace


def summary(dataset, **kwargs):
    """summarize heat_flux data in xarray DATASET

    Args:
        dataset(xarray.Dataset): the trajectory.dataset

    Returns:
        (pd.Dataframe, pd.Dataframe): One dataframe each for time resolved
              hf_acf, cumulative kappa, and frequency resolved spectra
    """

    assert keys.hf_acf in dataset
    assert keys.k_cum in dataset

    J1 = dataset[keys.hf_acf]
    Jw1 = get_fourier_transformed(J1).real

    js1 = Jw1.sum(axis=(1, 2))
    dct_freq = {keys.hf_power: js1}

    # get auxiliary heat_flux
    J2 = None
    if keys.hf_aux_acf in dataset:
        J2 = dataset[keys.hf_aux_acf]
        Jw2 = get_fourier_transformed(J2).real
        js2 = Jw2.sum(axis=(1, 2))
        dct_freq.update({keys.hf_aux_power: js2})

    # scalar
    kappa = dataset[keys.kappa_cumulative]
    kappa_scalar = xtrace(kappa) / 3
    J = xtrace(J1) / 3

    # time resolved
    d = {
        keys.hf_acf: J,
        keys.kappa_cumulative_scalar: kappa_scalar,
        keys.avalanche_function: dataset[keys.avalanche_function],
    }
    df_time = pd.DataFrame(d, index=dataset.time)

    # freq resolved
    df_freq = pd.DataFrame(dct_freq, index=Jw1.omega)

    return (df_time, df_freq)  #


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
