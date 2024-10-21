from argparse import ArgumentParser as argpars

import click
import numpy as np
import seaborn as sns
import xarray as xr
from matplotlib import pyplot as plt

from vibes import keys
from vibes.helpers.plotting import rc_params

plt.style.use(rc_params)


def plot_gk_summary(
    ds_gk: xr.Dataset,
    cmap: str = "colorblind",
    xlim: float = None,
    outfile: str = "greenkubo_summary.pdf",
):

    times_ps = ds_gk[keys.time] / 1000

    ks = ds_gk[keys.kappa]

    ks_flat = ks.stack(ab=("a", "b"))[::4].data
    k_mean = ks_flat.mean()
    k_err = (ks_flat.var() / (ks_flat.size)) ** 0.5

    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

    colors = sns.color_palette(cmap, n_colors=3)  # plt.get_cmap(cmap)

    cutoff_time = ds_gk[keys.time_cutoff] / 1000
    j_raw = ds_gk[keys.hf_acf]
    k_raw = ds_gk[keys.kappa_cumulative]
    j_filtered = ds_gk[keys.hf_acf_filtered]
    k_filtered = ds_gk[keys.kappa_cumulative_filtered]

    # diagonal values via stack
    k_raw_diag = k_raw.stack(ab=("a", "b"))[:, ::4]
    k_filtered_diag = k_filtered.stack(ab=("a", "b"))[:, ::4]

    k_total = k_raw_diag.mean(axis=1)
    k_total_filtered = k_filtered_diag.mean(axis=1)

    # Cartesian label
    labels = ["xx", "yy", "zz"]

    for ii in range(3):
        c = colors[ii]

        j = j_filtered[:, ii, ii]
        ax1.plot(times_ps, j, c=c, lw=2, label="$J_{"+f"{labels[ii]}""}$")

        # unfiltered kappa
        k = k_raw[:, ii, ii]
        ax2.plot(times_ps, k, c=c, alpha=0.5)

        k = k_filtered[:, ii, ii]
        ax2.plot(times_ps, k, c=c, label=r"$\kappa_{"+f"{labels[ii]}"+"}$")

        ta = cutoff_time[ii, ii]
        ax1.axvline(ta, c=c, lw=2)
        ax2.axvline(ta, c=c, lw=2)

        # unfiltered hfacf (and restore ylim)
        ylim = ax1.get_ylim()
        j = j_raw[:, ii, ii]
        ax1.plot(times_ps, j, c=c, lw=0.1, zorder=-1)
        ax1.set_ylim(ylim)

    # mean of k
    ax2.plot(times_ps, k_total, c="k", alpha=0.5)
    ax2.plot(times_ps, k_total_filtered, c="k", label="$\\kappa_{\\rm mean}$")

    ax1.axhline(0, c="k")
    ax1.set_ylim([j_filtered.min(), 1.2 * j_filtered.max()])

    # plot final kappa
    ax2.axhline(k_mean, c="k")
    ax2.fill_between(times_ps, k_mean + k_err, k_mean - k_err, color="k", alpha=0.1)

    ax1.set_ylabel("$J_{\\rm corr} (t)$")
    ax2.set_ylabel("$\\kappa (t)$ (W/mK)")
    ax2.set_xlabel("Time $t$ (ps)")

    kappa_str = f"$\\kappa$: {k_mean:.2f} +/- {k_err:.2f} W/mK"

    ax1.set_title(kappa_str)

    tmax = 3 * np.diag(cutoff_time).max()

    if xlim is None:
        xlim = tmax

    ax2.set_xlim([0, xlim])

    ax1.legend(loc="upper right", ncol=2)
    ax2.legend(loc="upper right", ncol=2)
    fig.tight_layout()

    if outfile is not None:
        fig.savefig(outfile)
        click.echo(f"..    green kubo summary plotted to {outfile}")

def main():
    """Main routine, deprecated since CLI"""
    parser = argpars(description="Plot summary for Green-Kubo formula")
    parser.add_argument("filename", help="greenkubo.nc")
    parser.add_argument("--outfile", default="greenkubo_summary.pdf")
    args = parser.parse_args()

    file = args.filename
    ds_gk = xr.load_dataset(file)
    plot_gk_summary(ds_gk, outfile=args.outfile)

if __name__ == "__main__":
    main()
