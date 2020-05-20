"""gather statistics about trajectory data"""
import numpy as np

# import pandas as pd
import xarray as xr
from ase.units import GPa

from vibes import keys
from vibes.helpers import talk, warn

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

    if len(series) < 1:
        return

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
    pprint(f"Temperature (last 1/2):  ", T_sum(series.iloc[len(series) // 2 :]))


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
    symbols = dataset.attrs["symbols"]
    usymbols = np.unique(symbols)

    # displacements
    dr = np.linalg.norm(dataset.displacements, axis=2)
    dr_mean = np.mean(dr, axis=1)
    dr_std = np.std(dr, axis=1)

    print()
    talk("Summarize Displacements", prefix="info")
    pprint(f"Avg. Displacement:", f"{dr.mean():.5} AA")
    pprint(f"Max. Displacement:", f"{dr.max():.5} AA")
    for sym in usymbols:
        mask = np.array(symbols) == sym
        # forces = dataset.forces[:, mask, :].data
        pprint(f"Avg. Displacement [{sym}]:", f"{dr[:, mask].mean():.5} AA")

    # forces
    forces = dataset.forces.data
    print()
    talk("Summarize Forces", prefix="info")
    pprint(f"Avg. Force:", f"{forces.mean():.5} eV/AA")
    pprint(f"Std. Force:", f"{forces.std():.5} eV/AA")

    for sym in usymbols:
        mask = np.array(symbols) == sym
        # forces = dataset.forces[:, mask, :].data
        pprint(f"Std. Force [{sym}]:", f"{forces[:, mask].std():.5} eV/AA")

    print()
    talk("Summarize Temperature", prefix="info")
    temperature(dataset.temperature)
    print()
    talk("Summarize Potential Energy", prefix="info")
    energy(dataset[keys.energy_potential])
    print()
    talk("Summarize Kinetic Pressure", prefix="info")
    pressure(dataset.pressure_kinetic)
    talk("Summarize Potential Pressure", prefix="info")
    pressure(dataset.pressure_potential)
    talk("Summarize Total Pressure (Kinetic + Potential)", prefix="info")
    pressure(dataset.pressure)

    # drift
    momenta = dataset.momenta.data
    momenta_time = np.sum(momenta, axis=1)
    momenta_mean = np.mean(abs(momenta_time), axis=0)
    rep = np.array2string(momenta_mean, precision=4)
    print()
    talk("Drift", prefix="info")
    pprint(f"Mean abs. Momentum:", f"{rep}")
    if any(momenta_mean > 1e-12):
        warn("Is is this drift problematic?", level=1)

    if plot:
        _keys = [
            keys.temperature,
            keys.energy_kinetic,
            keys.energy_potential,
            keys.pressure,
            keys.pressure_kinetic,
            keys.pressure_potential,
        ]
        df = dataset[_keys].to_dataframe()
        df.index /= 1000

        df["dr_mean"] = dr_mean
        df["dr_std"] = dr_std
        # df["dr"] = displacements
        # df["dr_std"] = dr_std

        plot_summary(df, **kwargs)
