"""compute and analyze heat fluxes"""
import xarray as xr

from ase import units
from hilde.helpers import Timer, warn


def get_velocities_data(trajectory, verbose=True):
    """extract velocties from TRAJECTORY  and return as xarray.DataArray

    Args:
        trajectory (Trajectory): list of atoms objects
    Returns:
        velocities (xarray.DataArray [N_t, N_a, 3])
    """
    timer = Timer("Get velocities from trajectory", verbose=verbose)

    metadata = {
        "time unit": "fs",
        "timestep": trajectory.timestep,
        "atoms": trajectory.symbols,
    }

    times = trajectory.times
    velocities = trajectory.velocities

    df = xr.DataArray(
        velocities,
        dims=["time", "atom", "coord"],
        coords={"time": times},
        name="velocities",
        attrs=metadata,
    )

    timer()

    return df


def get_pressure_data(trajectory, GPa=False, verbose=True):
    """extract pressure from TRAJECTORY  and return as xarray.DataArray

    Args:
        trajectory (Trajectory): list of atoms objects
    Returns:
        pressure (xarray.DataArray [N_t])
    """
    timer = Timer("Get pressure from trajectory", verbose=verbose)

    metadata = {
        "time unit": "fs",
        "timestep": trajectory.timestep,
        "atoms": trajectory.symbols,
        "GPa": GPa,
    }

    times = trajectory.times

    pressure = trajectory.pressure

    if GPa:
        pressure /= units.GPa

    df = xr.DataArray(
        pressure, dims=["time"], coords={"time": times}, name="pressure", attrs=metadata
    )

    timer()

    return df


def get_heat_flux_data(trajectory, return_avg=False):
    """compute heat fluxes from TRAJECTORY and return as xarray

    Args:
        trajectory: list of atoms objects WITH ATOMIC STRESS computed
        return_avg (bool): return flux per atom

    Returns:
        xarray.Dataset:
            heat_flux
            avg_heat_flux
    """
    traj_w_stresses = trajectory.with_stresses
    if not len(trajectory) == len(traj_w_stresses):
        warn("remove atoms from trajectory without stresses", level=1)
        trajectory = traj_w_stresses

    vec_index = ["time", "atom", "i"]
    scalar_index = "time"

    dataset = {
        "heat_flux": (vec_index, trajectory.heat_flux),
        "avg_heat_flux": (vec_index, trajectory.avg_heat_flux),
        "velocities": (vec_index, trajectory.velocities),
        "forces": (vec_index, trajectory.forces),
        "temperature": (scalar_index, trajectory.temperatures),
    }
    coords = {"time": trajectory.times}
    attrs = {
        "symbols": trajectory.symbols,
        "masses_dict": trajectory.masses_dict,
        "reference positions": trajectory.ref_positions,
    }

    ds = xr.Dataset(dataset, coords=coords, attrs=attrs)

    return ds
