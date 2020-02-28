"""`vibes info` backend"""
from pathlib import Path

from .misc import ClickAliasedGroup, click, complete_filenames


@click.command(cls=ClickAliasedGroup)
def info():
    """inform about content of a file"""


@info.command("geometry")
@click.argument("filename", default="geometry.in", type=complete_filenames)
@click.option("--format", default="aims", show_default=True)
@click.option("-t", "--symprec", default=1e-5, show_default=True)
@click.option("-v", "--verbose", is_flag=True)
@click.pass_obj
def geometry_info(obj, filename, format, symprec, verbose):
    """inform about a structure in a geometry input file"""
    from ase.io import read
    from vibes.structure.io import inform

    atoms = read(filename, format=format)

    verbosity = 1
    if verbose:
        verbosity = 2

    inform(atoms, symprec=symprec, verbosity=verbosity)


@info.command("md")
@click.argument("filename", default="trajectory.son", type=complete_filenames)
@click.option("-p", "--plot", is_flag=True, help="plot a summary")
@click.option("-w", "--write", is_flag=True, help="write Dataset to nc file")
@click.option("--avg", default=100, help="window size for running avg")
@click.option("-v", "--verbose", is_flag=True, help="be verbose")
def md_info(filename, plot, write, avg, verbose):
    """inform about content of a settings.in file"""
    import xarray as xr
    from .scripts.md_sum import md_sum
    from vibes.trajectory import analysis as al, reader

    file = Path(filename)

    if file.suffix in (".son", ".yaml", ".bz", ".gz"):
        trajectory = reader(filename)
        DS = trajectory.dataset
        if write:
            trajectory.write(file.parent / f"{file.stem}.nc")
    elif file.suffix in (".nc"):
        DS = xr.load_dataset(filename)
    elif file.suffix in (".log"):
        md_sum(filename, plot, avg, verbose)
    else:
        raise click.FileError(f"File format of {filename} not known.")

    click.echo(f"Dataset summary for {filename}:")
    al.summary(DS, plot=plot, avg=avg)


@info.command("phonopy")
@click.argument("filename", default="phonopy.in", type=complete_filenames)
@click.option("--write_supercell", is_flag=True, help="write the supercell to file")
@click.option("--format", default="aims", show_default=True)
def phonopy_info(filename, write_supercell, format):
    """inform about a phonopy calculation before it is started"""
    from .scripts.vibes_phonopy import preprocess

    preprocess(
        filename=None,
        settings_file=filename,
        dimension=None,
        format=format,
        write_supercell=write_supercell,
    )


@info.command("trajectory")
@click.argument("filename", default="trajectory.son", type=complete_filenames)
def trajectory_info(filename):
    """inform about content of trajectory file"""
    from vibes import son
    from vibes.settings import Settings

    metadata, _ = son.load(filename)

    click.echo(f"Summary of metadata in {filename}:\n")
    click.echo("Keys:")
    click.echo(f"  {list(metadata.keys())}\n")
    if "settings" in metadata:
        settings = Settings.from_dict(metadata["settings"])
        click.echo("Settings:")
        settings.print()


@info.command("netcdf")
@click.argument("file", type=complete_filenames)
def show_netcdf_file(file):
    """show contents of netCDF FILE"""
    import xarray as xr

    DS = xr.open_dataset(file)

    print(DS)


@info.command("hdf5")
@click.argument("file", type=complete_filenames)
@click.option("-v", "--verbose", is_flag=True)
def show_hdf5_file(file, verbose):
    """show contents of HDF5 FILE"""
    import pandas as pd

    click.echo(f"Summarize file {file}")
    with pd.HDFStore(file) as store:
        click.echo("Keys:")
        for k in store:
            click.echo(f"  {k}")

        if verbose:
            click.echo()
            for k in store:
                click.echo(f"Describe {k}")
                click.echo(store[k].describe())


@info.command(aliases=["gk"])
@click.argument("dataset", default="greenkubo.nc")
@click.option("-p", "--plot", is_flag=True, help="plot summary")
@click.option("--no_hann", is_flag=True)
@click.option("--logx", is_flag=True)
@click.option("--xlim", type=float, help="xlim range in ps")
@click.option("-avg", "--average", default=100, help="average window")
def greenkubo(dataset, plot, no_hann, logx, xlim, average):
    import xarray as xr
    from vibes import keys
    from vibes.green_kubo.analysis import summary, plot_summary

    DS = xr.load_dataset(dataset)

    (df_time, df_freq) = summary(DS)

    if plot:
        fig = plot_summary(
            df_time,
            df_freq,
            t_avalanche=DS.attrs[keys.time_avalanche],
            logx=logx,
            xlim=xlim,
            avg=average,
        )

        fname = Path(dataset).stem + "_summary.pdf"
        fig.savefig(fname, bbox_inches="tight")
        click.echo(f".. summary plotted to {fname}")


@info.command("vdos")
@click.argument("filename", default="trajectory.nc", type=complete_filenames)
@click.option("-o", "--output_filename", default="vdos.csv")
@click.option("-p", "--plot", is_flag=True, help="plot the DOS")
@click.option("--peak", type=float, help="height for peak detection", show_default=1)
@click.option("-mf", "--max_frequency", default=30.0, help="max. freq. in THz")
def velocity_autocorrelation(filename, output_filename, plot, peak, max_frequency):
    """write velocity autocorrelation function to output file"""
    import xarray as xr
    from vibes.green_kubo.velocities import get_vdos, simple_plot

    click.echo(f"Read {filename} and extract velocities")
    velocities = xr.open_dataset(filename).velocities

    vdos = get_vdos(velocities=velocities, hann=False, verbose=True)

    # sum atoms and coordinates
    df = vdos.real.sum(axis=(1, 2)).to_series()

    if plot:
        simple_plot(df, height=peak, max_frequency=max_frequency)

    click.echo(f".. write VDOS to {output_filename}")
    df.to_csv(output_filename, index_label="omega", header=True)
