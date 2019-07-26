"""compute and analyze heat fluxes"""
import xarray as xr

from ase import units
from hilde.helpers import warn
from . import Timer

time_index = "time"
vec_index = [time_index, "atom", "i"]
stress_index = [time_index, "i", "j"]


def _time_coords(trajectory):
    """return time as coords dict"""
    coords = {time_index: trajectory.times}
    return coords


def _metadata(trajectory, dct=None):
    """return metadata dictionary with defaults + custom dct"""
    attrs = {
        "time unit": "fs",
        "timestep": trajectory.timestep,
        "volume": trajectory.volume,
        "symbols": trajectory.symbols,
        "flattend reference positions": trajectory.ref_positions.flatten(),
    }

    if dct and isinstance(dct, dict):
        attrs.update(dct)

    return attrs


def get_positions_data(trajectory, verbose=True):
    """extract positions from TRAJECTORY  and return as xarray.DataArray

    Args:
        trajectory (Trajectory): list of atoms objects
    Returns:
        positions (xarray.DataArray [N_t, N_a, 3])
    """
    timer = Timer("Get positions from trajectory", verbose=verbose)

    df = xr.DataArray(
        trajectory.positions,
        dims=vec_index,
        coords=_time_coords(trajectory),
        name="positions",
        attrs=_metadata(trajectory),
    )

    timer()

    return df


def get_velocities_data(trajectory, verbose=True):
    """extract velocties from TRAJECTORY  and return as xarray.DataArray

    Args:
        trajectory (Trajectory): list of atoms objects
    Returns:
        velocities (xarray.DataArray [N_t, N_a, 3])
    """
    timer = Timer("Get velocities from trajectory", verbose=verbose)

    df = xr.DataArray(
        trajectory.velocities,
        dims=vec_index,
        coords=_time_coords(trajectory),
        name="velocities",
        attrs=_metadata(trajectory),
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

    unit = 1.0
    if GPa:
        unit = units.GPa

    extra_metadata = {"to_eV": unit}

    df = xr.DataArray(
        trajectory.pressure / unit,
        dims=[time_index],
        coords=_time_coords(trajectory),
        name="pressure",
        attrs=_metadata(trajectory, dct=extra_metadata),
    )

    timer()

    return df


def get_trajectory_data(trajectory):
    """Return trajectory data as xarray.Dataset

    Args:
        trajectory: list of atoms objects WITH ATOMIC STRESS computed
    Returns:
        xarray.Dataset:
            positions, velocities, forces, stress, pressure, temperature
    """

    # add velocities and pressure
    positions = get_positions_data(trajectory)
    velocities = get_velocities_data(trajectory)
    pressure = get_pressure_data(trajectory)

    dataset = {
        "positions": positions,
        "velocities": velocities,
        "forces": (vec_index, trajectory.forces),
        "stress": (stress_index, trajectory.stress),
        "pressure": pressure,
        "temperature": (time_index, trajectory.temperatures),
    }
    coords = _time_coords(trajectory)
    attrs = _metadata(trajectory)

    return xr.Dataset(dataset, coords=coords, attrs=attrs)


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

    # add velocities and pressure
    data = get_trajectory_data(trajectory)

    dataset = {
        "heat_flux": (vec_index, trajectory.heat_flux),
        "avg_heat_flux": (vec_index, trajectory.avg_heat_flux),
        "positions": data.positions,
        "velocities": data.velocities,
        "forces": data.forces,
        "pressure": data.pressure,
        "temperature": data.temperature,
    }
    coords = _time_coords(trajectory)
    attrs = _metadata(trajectory)

    return xr.Dataset(dataset, coords=coords, attrs=attrs)
