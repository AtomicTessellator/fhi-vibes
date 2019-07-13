"""hilde CLI utils"""

import scipy.signal as sl
import pandas as pd
from hilde.trajectory import reader
from .misc import click, AliasedGroup, complete_filenames


@click.command(cls=AliasedGroup)
def aiGK():
    """tools for (ab initio) Green Kubo, WORK IN PROGRESS"""


@aiGK.command(cls=AliasedGroup)
def autocorrelation():
    """utils for working with structures"""
    ...


@autocorrelation.command("velocity")
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