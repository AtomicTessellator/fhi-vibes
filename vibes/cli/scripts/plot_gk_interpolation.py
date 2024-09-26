from argparse import ArgumentParser as argpars
from pathlib import Path

import click
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from vibes import keys
from vibes.helpers.plotting import rc_params

plt.style.use(rc_params)


def plot_gk_interpolation(
    gk_ds: xr.Dataset,
    outfile: Path = "greenkubo_interpolation.pdf",
):
    """Plot summary for interpolation"""
    tau_sq = gk_ds[keys.mode_lifetime]

    # kappa = gk_ds[keys.kappa]
    kappa_ha = gk_ds[keys.kappa_ha]

    # correction = gk_ds[keys.interpolation_correction]

    # plot
    k_ai_r = gk_ds.heat_flux_autocorrelation_cumtrapz.stack(ab=("a", "b"))[:, ::4]
    _k_ha_q = gk_ds.heat_flux_harmonic_q_autocorrelation_cumtrapz
    k_ha_q = _k_ha_q.stack(ab=("a", "b"))[:, ::4]

    fig, ax = plt.subplots()
    ax.plot(k_ha_q.time, k_ai_r.mean(axis=1).data, label=r"$\kappa_{\rm aiGK}$")
    ax.plot(k_ha_q.time, k_ha_q.mean(axis=1).data, label=r"$\kappa_{\rm hm-q}$")
    km, kerr = k_ai_r.mean(axis=1), k_ai_r.std(axis=1) / 3 ** 0.5
    ax.fill_between(k_ai_r.time, km + kerr, km - kerr, alpha=0.25, color="C0")
    # k0 = np.diagonal(kappa).mean()
    k1 = np.diagonal(kappa_ha).mean()
    # k2 = np.diagonal(kappa).mean() + correction
    # k3 = np.diagonal(kappa).mean() * correction_factor
    # ax.axhline(k0, c="C0")
    # ax.axhline(k1, c="C1")
    # ax.axhline(k2, c="C0", ls="--")
    # ax.axhline(k3, c="C2", ls="--")
    ax.set_ylabel(r"$\kappa (t)$ (W/mK)")
    ax.set_xlabel("$t$ (fs)")

    kw = {"zorder": -1, "linestyle": "--"}
    cutoff_times = gk_ds.cutoff_time.stack(ab=("a", "b"))[::4]
    ax.axvline(cutoff_times.mean(), c="C0", **kw, label="Cutoff time")
    # for ct in cutoff_times:
    #     ax.axvline(ct, c="C0", lw=0.33, **kw)
    # ax.axvline(gk_ds.mode_lifetime.max(), c="C1", **kw)
    # ax.axvline(gk_ds.mode_lifetime_symmetrized.max(), c="C1", ls="--", **kw)

    ax.set_xlim(1e2)
    ax.set_xscale("log")

    ax.legend()
    fig.tight_layout()

    if outfile is not None:
        fig.savefig(outfile)
        click.echo(f".. interpolation summary plotted to {outfile}")

    # kappa harmonic
    fig, ax = plt.subplots()
    m = float(gk_ds.interpolation_fit_slope)
    nq = len(gk_ds.q_points) ** (1 / 3)
    ax.scatter(1 / nq, k1, marker="D", color="green", label=r"$\kappa_{\rm hm}$")

    # interpolation fit
    array = gk_ds[keys.interpolation_kappa_array]
    s = array.stack(sq=("a", "b"))[:, ::4].mean(axis=1).to_series()
    s.index = 1 / s.index
    s.plot(ax=ax, style="o", label=r"$\kappa_{\rm hm-int}$")

    # extrapolation
    y0 = float(gk_ds.interpolation_fit_intercept)
    x = np.linspace(0, 1.2 * 1 / nq, 3)
    ax.plot(x, x * m + y0, label="Extrapolation")
    ax.plot(0, y0, "X", clip_on=False, label=r"$\kappa_{\rm hm-bulk}$")
    ax.set_ylabel(r"$\kappa_{n_q}$ (W/mK)")
    ax.set_xlabel("$1/n_q$")
    ax.set_xlim(0)

    ax.legend()
    fig.tight_layout()

    if outfile is not None:
        _outfile = Path(outfile).stem + "_fit" + Path(outfile).suffix
        fig.savefig(_outfile)
        click.echo(f".. interpolation summary plotted to {_outfile}")

    # lifetimes
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    kw = {"color": "black", "alpha": 0.05}

    x = gk_ds.time.data / 1000  # in ps

    y1 = []
    y2 = []
    for sq in np.ndindex(tau_sq.shape):
        tau = tau_sq[sq]
        if np.isnan(tau):
            continue
        y1.append(gk_ds.mode_energy_autocorrelation.data[:, sq[0], sq[1]])
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

def main():
    """Main routine, deprecated since CLI"""
    parser = argpars(description="Read geometry and use spglib to refine")
    parser.add_argument("filename", help="greenkubo.nc")
    parser.add_argument("--outfile", default="greenkubo_interpolation.pdf")
    args = parser.parse_args()

    file = args.filename
    ds_gk = xr.load_dataset(file)
    plot_gk_interpolation(ds_gk, outfile=args.outfile)

if __name__ == "__main__":
    main()
