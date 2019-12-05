"""test VDOS functionality"""

from pathlib import Path
import numpy as np
import scipy.signal as sl

from ase.io import read
from vibes.trajectory import reader
from vibes.green_kubo.velocities import get_vdos
from vibes.tdep.wrapper import parse_tdep_forceconstant
from vibes.harmonic_analysis import dynamical_matrix as dm

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


def test_vdos():
    traj = reader(parent / "trajectory.son.bz2")

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
