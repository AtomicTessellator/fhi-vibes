"""utilities for working with trajectories"""
import numpy as np


def clean_pressure(series):
    """remove 0 pressure entries from a SERIES of pressures

    Reference:
        https://pandas.pydata.org/pandas-docs/stable/user_guide/missing_data.html
    """
    series[abs(series) < 1e-20] = np.nan
    return series
