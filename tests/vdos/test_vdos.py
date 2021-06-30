"""test VDOS functionality"""

from pathlib import Path

import numpy as np
import scipy.signal as sl
from ase.io import read

from vibes.dynamical_matrix import get_frequencies
from vibes.green_kubo.velocities import get_vdos
from vibes.io import parse_force_constants
from vibes.trajectory import reader


parent = Path(__file__).parent


def test_parse_force_constants():
    # frequencies from force constants
    fc = parse_force_constants(
        parent / "FORCE_CONSTANTS_tdep",
        primitive=parent / "geometry.in.primitive",
        supercell=parent / "geometry.in.supercell",
        two_dim=True,
    )

    return fc


def test_frequencies_from_force_constants():
    fc = test_parse_force_constants()
    sc = read(parent / "geometry.in.supercell", format="aims")

    freqs = get_frequencies(fc, masses=sc.get_masses())

    return freqs


def test_vdos(traj_file="trajectory.son", vdos_file="v.nc", ref_file="ref_vdos.csv"):
    traj = reader(parent / traj_file)

    df_vdos = get_vdos(traj.dataset.velocities, npad=0).real

    # get analytical frequencies
    freqs = test_frequencies_from_force_constants()

    # convert to pandas.Series
    ds = df_vdos.sum(axis=(1, 2)).to_series()

    # compare peak positions and analytical frequencies
    peaks = ds.iloc[sl.find_peaks(ds.to_numpy())[0]].index

    unique_freqs = np.unique(np.round(freqs.real, decimals=3))[1:]

    for peak, freq in zip(peaks, unique_freqs):
        assert abs(peak - freq) / peak < 0.1, (peak, freq)

    # check zero padding
    df_vdos_pad = get_vdos(traj.dataset.velocities, npad=100000).real
    ds_pad = df_vdos_pad.sum(axis=(1, 2)).to_series()

    peaks_pad = ds_pad.iloc[sl.find_peaks(ds_pad, height=0.2)[0]].index

    for peak, freq in zip(peaks_pad, unique_freqs):
        diff = abs(peak - freq) / peak
        assert diff < 0.1, (peak, freq)


if __name__ == "__main__":
    test_frequencies_from_force_constants()
    test_parse_force_constants()
    test_vdos()
