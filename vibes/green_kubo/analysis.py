"""gather statistics about trajectory data"""
import pandas as pd

from vibes import keys
from vibes.correlation import get_autocorrelation_exponential
from vibes.fourier import get_fourier_transformed
from vibes.helpers.xarray import xtrace


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


def plot_summary(df_time, df_freq, t_avalanche=None, avg=50, logx=True, xlim=None):
    """plot a summary of the data in DATAFRAME"""
    import matplotlib
    from matplotlib import pyplot as plt
    from vibes.helpers.plotting import tableau_colors as tc

    matplotlib.use("pdf")

    try:
        import seaborn as sns

        sns.set_style("whitegrid")
        sns.set_palette("colorblind")
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
    kw_avalanche = {"color": tc[0], "linestyle": "--", "linewidth": 2}
    avg_kw = {"linewidth": 3, "color": "k"}
    fig_kw = {
        "figsize": (11.69, 8.27),
        "gridspec_kw": {"width_ratios": [3, 2]},
        "sharex": "col",
    }
    kw_roll = {"window": avg, "min_periods": 0, "center": True}
    # kw_exp1 = {"min_periods": 0, "center": False}
    # kw_exp2 = {"min_periods": 0}  # "center": True}

    df_time.index /= 1000

    jc = df_time[keys.hf_acf] / df_time[keys.hf_acf].iloc[0]
    kc = df_time[keys.kappa_cumulative_scalar]
    js = df_freq[keys.hf_power]

    # estimate correlation time
    e = get_autocorrelation_exponential(jc, ps=True, verbose=False)

    fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(nrows=2, ncols=2, **fig_kw)

    jc.rolling(5, min_periods=0).mean().plot(ax=ax11)

    e.plot(ax=ax11, zorder=10)

    # plot
    # HFACF
    jc.plot(ax=ax11, **plot_kw)
    jc.rolling(**kw_roll).mean().plot(ax=ax11, **avg_kw)

    # avalanche
    twin = ax11.twinx()
    fig.align_ylabels()
    # twin.set_yticks([]), twin.set_yticklabels([])
    f = df_time[keys.avalanche_function]
    f.plot(ax=twin)
    twin.grid(False)
    if t_avalanche is not None:
        ax11.axvline(t_avalanche / 1000, **kw_avalanche)

    # Kappa
    kc.plot(ax=ax21, **{**avg_kw, "color": tc[1], "linewidth": 2})
    ax21.axvline(t_avalanche, **kw_avalanche)

    # plot spectra
    js.plot(ax=ax12, **plot_kw)
    js.rolling(**kw_roll).mean().plot(ax=ax12, **avg_kw)

    if xlim:
        ax21.set_xlim((-1, xlim))
    else:
        ax21.set_xlim((-1, kc.index.max()))
    ax21.set_ylim((0, 1.1 * kc.max()))
    if logx:
        ax21.set_xscale("log")
        ax21.set_xlim(kc.index[2], kc.index.max())

    # aux flux
    if keys.hf_aux_power in df_freq:
        js_aux = df_freq[keys.hf_aux_power]
        js_aux.plot(ax=ax22, **plot_kw)
        js_aux.rolling(**kw_roll).mean().plot(ax=ax22, **avg_kw)
        linthreshy = 0.1 * max(1, js_aux.max())
        ax22.set_yscale("symlog", linthreshy=linthreshy)
        ax22.axhline(linthreshy, ls="--", color="purple")
        ax22.set_xlim((-1, js_aux.index.max()))
        ax22.set_xlabel("Omega (THz)")
        ax22.set_title(r"$\vert J_\mathrm{aux} (\omega)\vert^2$")
    else:
        ax22.set_title(r"$J_\mathrm{aux} (\omega)$ missing.")

    # titles and labels
    fig.suptitle(f"Heat Flux Overview (window: {avg})")
    ax11.set_title(r"$\langle J (t) J(0) \rangle ~/~ \langle J (0) J(0) \rangle$")
    ax21.set_title(r"$\kappa (t)$ (W/mK)")
    ax21.set_xlabel("Time $t$ (ps)")
    ax12.set_title(r"$\vert J (\omega)\vert^2$")
    return fig
