"""hilde CLI utils"""

import scipy.signal as sl
import xarray as xr
from hilde.trajectory import reader
from .misc import click, AliasedGroup, complete_filenames


@click.command(cls=AliasedGroup)
def aiGK():
    """tools for (ab initio) Green Kubo, WORK IN PROGRESS"""


@aiGK.command(cls=AliasedGroup)
def autocorrelation():
    """utils for working with structures"""
    ...


@autocorrelation.command("vdos")
@click.argument("filename", default="velocities.nc", type=complete_filenames)
@click.option("-o", "--output_filename", default="vdos.csv")
@click.option("-p", "--plot", is_flag=True, help="plot the DOS")
def velocity_autocorrelation(filename, output_filename, plot):
    """write velocity autocorrelation function to output file"""
    from hilde.green_kubo.velocities import get_vdos, simple_plot

    velocities = xr.open_dataarray(filename)

    vdos = get_vdos(velocities=velocities, verbose=True)

    # sum atoms and coordinates
    df = vdos.real.sum(axis=(1, 2)).to_series()

    click.echo(f".. write VDOS to {output_filename}")
    df.to_csv(output_filename, index_label="omega", header=True)

    if plot:
        simple_plot(df)


@click.argument("filename", type=complete_filenames)
@click.option("-o", "--output_filename", default="velocities.csv")
def velocity_autocorrelation(filename, output_filename):
    """write velocity autocorrelation function to output file"""

    traj = reader(filename)

    times = traj.times
    e_kin = []
    velocities = []
    for atoms in traj:
        v = atoms.get_velocities()
        velocities.append(v)
        e = atoms.get_kinetic_energy()
        e_kin.append(e)

    assert len(times) == len(velocities)
    assert len(times) == len(e_kin)

    df = pd.DataFrame({"e_kin": e_kin}, index=times)
    C_e = sl.correlate(e_kin, e_kin)[len(e_kin) - 1 :] / len(e_kin)
    df["e_kin_corr"] = C_e

    # vv(\tau) = \sum_i v_i (tau) v_i (0)

    df.to_csv(output_filename, index_label="time")