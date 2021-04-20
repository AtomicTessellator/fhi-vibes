import click
import numpy as np
import seaborn as sns
import xarray as xr
from matplotlib import pyplot as plt
from vibes import keys
from vibes.helpers.plotting import rc_params


plt.style.use(rc_params)


def main(
    file: str = "greenkubo.nc",
    cmap: str = "colorblind",
    xlim: float = None,
    outfile: str = "gk_summary.pdf",
):
    if isinstance(file, xr.Dataset):
        ds_gk = file
    else:
        ds_gk = xr.load_dataset(file)

    ds_gk[keys.time] = ds_gk[keys.time] / 1000

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

    for ii in range(3):
        c = colors[ii]

        j = j_filtered[:, ii, ii]
        j.to_series().plot(ax=ax1, c=c, lw=2)

        # unfiltered kappa
        k = k_raw[:, ii, ii]
        k.to_series().plot(ax=ax2, c=c, alpha=0.75)

        k = k_filtered[:, ii, ii]
        k.to_series().plot(ax=ax2, c=c)

        ta = cutoff_time[ii, ii]
        ax1.axvline(ta, c=c, lw=2)
        ax2.axvline(ta, c=c, lw=2)

        # unfiltered hfacf (and restore ylim)
        ylim = ax1.get_ylim()
        j = j_raw[:, ii, ii]
        j.to_series().plot(ax=ax1, c=c, lw=0.1, zorder=-1)
        ax1.set_ylim(ylim)

    # mean of k
    k_total.to_series().plot(ax=ax2, c="k", alpha=0.75)
    k_total_filtered.to_series().plot(ax=ax2, c="k")

    ax1.axhline(0, c="k")
    ax1.set_ylim([j_filtered.min(), 1.2 * j_filtered.max()])

    # plot final kappa
    ax2.axhline(k_mean, c="k")
    ax2.fill_between(j.time, k_mean + k_err, k_mean - k_err, color="k", alpha=0.1)

    ax1.set_ylabel("$J_{\\rm corr} (t)$")
    ax2.set_ylabel("$\\kappa (t)$")
    ax2.set_xlabel("Time $t$ (ps)")

    kappa_str = f"$\\kappa$: {k_mean:.2f} +/- {k_err:.2f} W/mK"

    ax1.set_title(kappa_str)

    tmax = 3 * np.diag(cutoff_time).max()

    if xlim is None:
        xlim = tmax

    ax2.set_xlim([0, xlim])

    if outfile is not None:
        fig.savefig(outfile)
        click.echo(f"..    green kubo summary plotted to {outfile}")


if __name__ == "__main__":
    import typer

    typer.run(main)
