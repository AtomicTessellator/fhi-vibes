from pathlib import Path

import click
import matplotlib.transforms as transforms
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from vibes.helpers.plotting import rc_params


plt.style.use(rc_params)


def main(file: str = "greenkubo.nc", outfile: Path = "gk_summary_interpolation.png"):
    """plot summary for interpolation"""
    if isinstance(file, xr.Dataset):
        DS = file
    else:
        DS = xr.load_dataset(file)

    tau_sq = DS.mode_lifetime

    kappa = DS.thermal_conductivity
    kappa_ha = DS.thermal_conductivity_harmonic

    correction = DS.interpolation_correction

    tau_sq.stack(sq=("s", "q")).to_series().describe()

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

    fig.savefig(outfile)
    click.echo(f".. interpolation summary plotted to {outfile}")

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
        y1.append(DS.mode_occupation_acf.data[:, sq[0], sq[1]])
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
    ax1.set_ylabel(r"$g(t)$", rotation=0)

    fig.suptitle("Mode occupation autocorrelation $g(t)$")
    ax1.set_title("Simulation")
    ax2.set_title("Analytic")

    outfile = Path(outfile).stem + "_lifetimes" + Path(outfile).suffix
    fig.savefig(outfile)
    click.echo(f"..      lifetime summary plotted to {outfile}")


if __name__ == "__main__":
    import typer

    typer.run(main)
