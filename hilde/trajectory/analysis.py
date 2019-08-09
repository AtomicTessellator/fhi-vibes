"""gather statistics about trajectory data"""
# import numpy as np
# import pandas as pd
import xarray as xr
from ase.units import GPa
from hilde.helpers import talk
from .plotting import plot_summary


def np_thirds(array):
    """subdivide numpy or xarray array into thirds"""
    len_3 = len(array) // 3
    sub_1 = array[:len_3]
    sub_2 = array[len_3 : 2 * len_3]
    sub_3 = array[2 * len_3 :]

    return sub_1, sub_2, sub_3


def pd_thirds(series):
    """subdivide pandas.Series/DataFrame into thirds"""
    len_3 = len(series) // 3
    sub_1 = series.iloc[:len_3]
    sub_2 = series.iloc[len_3 : 2 * len_3]
    sub_3 = series.iloc[2 * len_3 :]

    return sub_1, sub_2, sub_3


def pprint(msg1, msg2, width1=25):
    """pretty print with fixed field size for first message"""
    print(f"{msg1:{width1}s} {msg2}")


def e_sum(e):
    return f"{e.mean():12.3f} +/- {e.std():12.4f} eV"


def T_sum(T):
    return f"{T.mean():12.3f} +/- {T.std():12.4f} K"


def p_sum(p, to_GPa=True):
    """print summary for pressure (slice)"""
    unit = "eV/AA**3"
    p = p.copy()
    if to_GPa:
        p /= GPa
        unit = "GPa"
    return f"{p.mean():12.6f} +/- {p.std():10.6f} {unit}"


def pressure(series, interpolate=None):
    """summarize pressure from MD

    Args:
        series: Series/Dataframe representing pressure
        interpolate: use `interpolate` to deal with missing values
    """

    if isinstance(series, xr.core.dataarray.DataArray):
        series = series.to_series()

    # remove zeros
    len_orig = len(series)
    if interpolate:
        series = series.interpolate(interpolate)
    else:
        series = series.dropna()

    # time in ps
    time = series.index / 1000

    # thirds
    # thirds = pd_thirds(series)

    msg = f"{time[-1] - time[0]:8.3f} ps ({len(time)} of {len_orig} steps)"
    pprint("Simulation time:", msg)
    pprint("Pressure:", p_sum(series))
    pprint("Pressure (last 1/2):", p_sum(series.iloc[len(series) // 2 :]))
    pprint("Pressure (last 1/2):", p_sum(series.iloc[len(series) // 2 :], to_GPa=False))
    # pprint(f"Pressure ({ii+1}st 1/3):", p_sum(p, to_GPa=False))


def temperature(series):
    """summarize temperature from MD"""

    if isinstance(series, xr.core.dataarray.DataArray):
        series = series.to_series()

    # time in ps
    time = series.index / 1000

    # thirds
    thirds = pd_thirds(series)

    msg = f"{time[-1] - time[0]:8.3f} ps ({len(time)} steps)"
    pprint("Simulation time:", msg)
    pprint("Temperature:", T_sum(series))
    for ii, T in enumerate(thirds):
        pprint(f"Temperature ({ii+1}st 1/3):", T_sum(T))


def energy(series):
    """summarize energies from MD"""

    if isinstance(series, xr.core.dataarray.DataArray):
        series = series.to_series()

    # time in ps
    time = series.index / 1000

    # thirds
    thirds = pd_thirds(series)

    msg = f"{time[-1] - time[0]:8.3f} ps ({len(time)} steps)"
    pprint("Simulation time:", msg)
    pprint("Pot. Energy:", e_sum(series))
    for ii, e in enumerate(thirds):
        pprint(f"Pot. Energy ({ii+1}st 1/3):", e_sum(e))


def summary(dataset, plot=False, **kwargs):
    """summarize MD data in xarray DATASET"""

    print()
    talk("Summarize Temperature", prefix="info")
    temperature(dataset.temperature)
    print()
    talk("Summarize Potential Energy", prefix="info")
    energy(dataset.potential_energy)
    print()
    talk("Summarize Pressure", prefix="info")
    pressure(dataset.pressure)

    if plot:
        df = dataset[
            ["temperature", "kinetic_energy", "potential_energy", "pressure"]
        ].to_dataframe()
        df.index /= 1000
        plot_summary(df, **kwargs)