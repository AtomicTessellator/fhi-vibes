"""vibes CLI utils"""

# import scipy.signal as sl
# import xarray as xr
# from vibes.trajectory import reader
from .misc import click, AliasedGroup, complete_filenames


@click.command(cls=AliasedGroup)
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


# @click.argument("filename", type=complete_filenames)
# @click.option("-o", "--output_filename", default="velocities.csv")
# def velocity_autocorrelation(filename, output_filename):
#     """write velocity autocorrelation function to output file"""
#
#     traj = reader(filename)
#
#     times = traj.times
#     e_kin = []
#     velocities = []
#     for atoms in traj:
#         v = atoms.get_velocities()
#         velocities.append(v)
#         e = atoms.get_kinetic_energy()
#         e_kin.append(e)
#
#     assert len(times) == len(velocities)
#     assert len(times) == len(e_kin)
#
#     df = pd.DataFrame({"e_kin": e_kin}, index=times)
#     C_e = sl.correlate(e_kin, e_kin)[len(e_kin) - 1 :] / len(e_kin)
#     df["e_kin_corr"] = C_e
#
#     # vv(\tau) = \sum_i v_i (tau) v_i (0)
#
#     df.to_csv(output_filename, index_label="time")
