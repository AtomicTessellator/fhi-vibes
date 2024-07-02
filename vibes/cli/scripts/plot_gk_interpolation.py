from pathlib import Path

import click
import matplotlib.transforms as transforms
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from vibes import keys
from vibes.helpers.plotting import rc_params


plt.style.use(rc_params)


def main(
    file: str = "greenkubo.nc", outfile: Path = "greenkubo_summary_interpolation.pdf"
):
    """plot summary for interpolation"""
    if isinstance(file, xr.Dataset):
        DS = file
    else:
        DS = xr.load_dataset(file)

    tau_sq = DS[keys.mode_lifetime]

    kappa = DS[keys.thermal_conductivity]
    kappa_ha = DS[keys.thermal_conductivity_harmonic]

    correction = DS[keys.interpolation_correction]

    # plot
    k_ai_r = DS.heat_flux_acf_integral.stack(ab=("a", "b"))[:, ::4]
    k_ha_q = DS.heat_flux_harmonic_q_acf_integral.stack(ab=("a", "b"))[:, ::4]

    fig, ax = plt.subplots()
    ax = k_ai_r.mean(axis=1).to_series().plot()
    km, kerr = k_ai_r.mean(axis=1), k_ai_r.std(axis=1) / 3 ** 0.5
    ax.fill_between(k_ai_r.time, km + kerr, km - kerr, alpha=0.25, color="C0")
    k_ha_q.mean(axis=1).to_series().plot(ax=ax)
    k0 = np.diagonal(kappa).mean()
    k1 = np.diagonal(kappa_ha).mean()
    k2 = np.diagonal(kappa).mean() + correction
    # k3 = np.diagonal(kappa).mean() * correction_factor
    ax.axhline(k0, c="C0")
    ax.axhline(k1, c="C1")
    ax.axhline(k2, c="C0", ls="--")
    # ax.axhline(k3, c="C2", ls="--")
    ax.set_ylabel(r"$\kappa (t)$ (W/mK)")
    ax.set_xlabel("$t$ (fs)")

    x = 0.05
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    kw = {"transform": trans, "va": "top"}
    ax.text(x, 0.99 * k0, "ai", color="C0", **kw)
    ax.text(x, 0.99 * k1, "ha-q", color="C1", **kw)
    ax.text(x, 0.99 * k2, "ai-ext.", color="C0", **kw)

    kw = {"zorder": -1}
    cutoff_times = DS.cutoff_time.stack(ab=("a", "b"))[::4]
    ax.axvline(cutoff_times.mean(), c="C0", **kw)
    for ct in cutoff_times:
        ax.axvline(ct, c="C0", lw=0.33, **kw)
    ax.axvline(DS.mode_lifetime.max(), c="C1", **kw)
    # ax.axvline(DS.mode_lifetime_symmetrized.max(), c="C1", ls="--", **kw)

    ax.set_xlim([1e2, 1e5])
    ax.set_xscale("log")

    fig = ax.get_figure()

    if outfile is not None:
        fig.savefig(outfile)
        click.echo(f".. interpolation summary plotted to {outfile}")

    # interpolation fit
    fig, ax = plt.subplots()
    array = DS[keys.interpolation_kappa_array]
    s = array.stack(sq=("a", "b"))[:, ::4].mean(axis=1).to_series()
    s.index = 1 / s.index
    s.plot(ax=ax, style=".-")

    nq = len(DS.q_points) ** (1 / 3)

    m = float(DS.interpolation_fit_slope)
    ax.scatter(1 / nq, k1, marker="D", color="green")
    x = np.linspace(0, 1.2 * 1 / nq, 3)
    ax.plot(x, (x - 1 / nq) * m + k1, zorder=-1)
    ax.set_ylabel(r"$\kappa_{n_q}$ (W/mK)")
    ax.set_xlabel("$1/n_q$")
    ax.set_xlim([0, 0.5])

    if outfile is not None:
        _outfile = Path(outfile).stem + "_fit" + Path(outfile).suffix
        fig.savefig(_outfile)
        click.echo(f".. interpolation summary plotted to {_outfile}")

    # lifetimes
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    kw = {"color": "black", "alpha": 0.05}

    x = DS.time.data / 1000  # in ps

    y1 = []
    y2 = []
    for sq in np.ndindex(tau_sq.shape):
        tau = tau_sq[sq]
        if np.isnan(tau):
            continue
        y1.append(DS.mode_energy_acf.data[:, sq[0], sq[1]])
        y2.append(np.exp(-x * 1000 / float(tau_sq[sq])))

    # plot in segments (memory)
    shape = (*np.shape(y1), 2)
    segments1 = np.zeros(shape)
    segments1[:, :, 0] = x
    segments1[:, :, 1] = y1
    segments2 = segments1.copy()
    segments2[:, :, 1] = y2

    lc1 = LineCollection(segments1, **kw)
    lc2 = LineCollection(segments2, **kw)
    ax1.add_collection(lc1)
    ax2.add_collection(lc2)

    ylim = [0.09, 1]
    yticks = [0.1, 1]
    for ax in (ax1, ax2):
        ax.set_xlim([0, 5])
        ax.set_xlabel("Time (ps)")
    ax1.set_ylim(ylim)
    ax1.set_yscale("log")
    ax1.set_yticks(yticks)
    ax1.set_yticks(np.arange(0.1, 1, 0.1), minor=True)
    ax1.set_yticklabels(yticks)
    ax1.set_yticklabels([], minor=True)
    ax1.set_ylabel(r"$G_s(t)$", rotation=0)

    fig.suptitle("Mode energy autocorrelation")
    ax1.set_title("Simulation")
    ax2.set_title("Analytic")

    if outfile is not None:
        _outfile = Path(outfile).stem + "_lifetimes" + Path(outfile).suffix
        fig.savefig(_outfile)
        click.echo(f"..      lifetime summary plotted to {_outfile}")

    # mode resolved
    # plot this component
    ab = [0,2]
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, sharex=True, figsize=(10, 5))
    kw = {"color": "black", "alpha": 0.5}

    x = DS.time.data / 1000  # in ps

    y1 = [] # diagonal
    y2 = [] # sum over all off-diagonal
    y3 = [] # kappa diagonal
    y4 = [] # kappa off-diag
    y1_lim_max = 0
    y3_lim_max = 0
    for s in np.arange(DS.heat_flux_harmonic_q_tsa_acf.shape[1]):
        y1.append(DS.heat_flux_harmonic_q_tsa_acf.data[:, s, ab[0], s, ab[1]])
        _y1_max = np.abs(DS.heat_flux_harmonic_q_tsa_acf.data[:, s, ab[0], s, ab[1]]).max()
        if _y1_max > y1_lim_max:
            y1_lim_max = _y1_max
        hf_ha_q_offdiag = DS.heat_flux_harmonic_q_tsa_acf.sum(axis=3)[:, s, ab[0], ab[1]] - \
                                DS.heat_flux_harmonic_q_tsa_acf[:, s, ab[0], s, ab[1]]
        y2.append(hf_ha_q_offdiag)
        y3.append(DS.heat_flux_harmonic_q_tsa_acf_integral.data[:, s, ab[0], s, ab[1]])
        _y3_max = np.abs(DS.heat_flux_harmonic_q_tsa_acf_integral.data[:, s, ab[0], s, ab[1]]).max()
        if _y3_max > y3_lim_max:
            y3_lim_max = _y3_max
        hf_ha_q_integral_offdiag = DS.heat_flux_harmonic_q_tsa_acf_integral.sum(axis=3)[:, s, ab[0], ab[1]] - \
                                DS.heat_flux_harmonic_q_tsa_acf_integral[:, s, ab[0], s, ab[1]]
        y4.append(hf_ha_q_integral_offdiag)
    
    # plot in segments (memory)
    shape = (*np.shape(y1), 2)
    segments1 = np.zeros(shape)
    segments1[:, :, 0] = x
    segments1[:, :, 1] = y1
    segments2 = segments1.copy()
    segments2[:, :, 1] = y2
    segments3 = segments1.copy()
    segments3[:, :, 1] = y3
    segments4 = segments1.copy()
    segments4[:, :, 1] = y4

    lc1 = LineCollection(segments1, **kw)
    lc2 = LineCollection(segments2, **kw)
    lc3 = LineCollection(segments3, **kw)
    lc4 = LineCollection(segments4, **kw)
    ax1.add_collection(lc1)
    ax2.add_collection(lc2)
    ax3.add_collection(lc3)
    ax4.add_collection(lc4)

    y1lim = [y1_lim_max*-1.2, y1_lim_max*1.2]
    y3lim = [y3_lim_max*-1.2, y3_lim_max*1.2]
    #yticks = [-0.005, 0.005]
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlim([0.01, x.max() * 1.2])
        ax.set_xscale("log")
        ax.set_xlabel("Time (ps)")
    ax1.set_ylim(y1lim)
    ax2.set_ylim(y1lim)
    ax3.set_ylim(y3lim)
    ax4.set_ylim(y3lim)
    #ax1.set_yscale("log")
    #ax1.set_yticks(yticks)
    #ax1.set_yticks(np.arange(0.1, 1, 0.1), minor=True)
    #ax1.set_yticklabels(yticks)
    #ax1.set_yticklabels([], minor=True)
    ax1.set_ylabel(r"$G_s(t)$", rotation=0)
    ax3.set_ylabel(r"$kappa_s(t)$", rotation=0)

    fig.suptitle("Harmonic mode-resolve correlation")
    ax1.set_title("Diagonal")
    ax2.set_title("Off-diagonal")

    if outfile is not None:
        _outfile = Path(outfile).stem + f"_ha_q_{ab[0]}{ab[1]}" + Path(outfile).suffix
        fig.savefig(_outfile)
        click.echo(f".. harmonic mode-resolve plotted to {_outfile}")


if __name__ == "__main__":
    import typer

    typer.run(main)
