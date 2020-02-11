"""vibes CLI utils"""
from .misc import ClickAliasedGroup, click, complete_filenames


@click.command(cls=ClickAliasedGroup)
def aiGK():
    """tools for (ab initio) Green Kubo, WORK IN PROGRESS"""


# @aiGK.command(cls=AliasedGroup)
# def autocorrelation():
#     """utils for working with structures"""
#     ...


@aiGK.command("vdos")
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

    vdos = get_vdos(velocities=velocities, verbose=True)

    # sum atoms and coordinates
    df = vdos.real.sum(axis=(1, 2)).to_series()

    click.echo(f".. write VDOS to {output_filename}")
    df.to_csv(output_filename, index_label="omega", header=True)

    if plot:
        simple_plot(df, height=peak, max_frequency=max_frequency)
