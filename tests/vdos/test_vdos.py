"""test VDOS functionality"""

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.signal as sl
import xarray as xr
from ase.io import read
from matplotlib import pyplot as plt

from vibes.green_kubo.velocities import get_vdos, get_velocity_autocorrelation
from vibes.harmonic_analysis import dynamical_matrix as dm
from vibes.tdep.wrapper import parse_tdep_forceconstant
from vibes.trajectory import reader

parent = Path(__file__).parent


def test_parse_force_constants():
    # frequencies from force constants
    fc = parse_tdep_forceconstant(
        parent / "infile.forceconstant",
        parent / "geometry.in.primitive",
        parent / "geometry.in.supercell",
        two_dim=True,
        format="aims",
    )

    return fc


def test_frequencies_from_force_constants():
    fc = test_parse_force_constants()
    sc = read(parent / "geometry.in.supercell", format="aims")

    freqs = dm.get_frequencies(fc, masses=sc.get_masses())

    return freqs


def test_vdos(
    traj_file="trajectory.son.bz2", vdos_file="v.nc", ref_file="ref_vdos.csv"
):
    traj = reader(parent / traj_file)

    df_vdos = get_vdos(trajectory=traj)

    # get analytical frequencies
    freqs = test_frequencies_from_force_constants()

    # convert to pandas.Series
    ds = df_vdos.sum(axis=(1, 2)).to_series()

    # compare peak positions and analytical frequencies
    peaks = ds.iloc[sl.find_peaks(ds.to_numpy().real)[0]].index

    unique_freqs = np.unique(np.round(freqs.real, decimals=3))[1:]

    for peak, freq in zip(peaks, unique_freqs):
        assert abs(peak - freq) / peak < 0.1, (peak, freq)

    # compare to ref
    velocities = xr.load_dataarray(parent / vdos_file)

    vdos = get_vdos(velocities=velocities).sum(axis=(1, 2)).to_series().abs()

    vdos_ref = pd.read_csv(parent / ref_file, index_col=vdos.index.name, squeeze=True)

    assert (vdos - vdos_ref).std() < 1e-12


if __name__ == "__main__":
    test_frequencies_from_force_constants()
    test_parse_force_constants()
    test_vdos()
