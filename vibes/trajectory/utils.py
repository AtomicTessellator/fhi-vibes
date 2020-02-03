"""utilities for working with trajectories"""
import numpy as np

from vibes.helpers import Timer
from vibes.helpers import talk as _talk


_prefix = "trajectory"

Timer.prefix = _prefix


def talk(msg):
    """wrapper for `utils.talk` with prefix"""
    return _talk(msg, prefix=_prefix)


def get_hashes_from_trajectory_file(trajectory, verbose=False):
    """return all hashes from trajectory"""
    from .io import reader

    try:
        traj = reader(trajectory, verbose=verbose)
    except (FileNotFoundError, KeyError):
        return []

    return traj.get_hashes()


def clean_pressure(series):
    """remove 0 pressure entries from a SERIES of pressures

    Reference:
        https://pandas.pydata.org/pandas-docs/stable/user_guide/missing_data.html
    """
    series[abs(series) < 1e-20] = np.nan
    return series
