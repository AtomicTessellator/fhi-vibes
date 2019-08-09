"""compute and analyze heat fluxes"""
import xarray as xr

from ase import units
from hilde.helpers import warn
from hilde.structure.misc import get_sysname
from .utils import clean_pressure
from . import Timer

time_dims = "time"
vec_dims = [time_dims, "I", "a"]
stress_dims = [time_dims, "a", "b"]
kappa_dims = [time_dims, "I", "J", "a", "b"]


def _time_coords(trajectory):
    """return time as coords dict"""
    coords = {time_dims: trajectory.times}
    return coords


def _metadata(trajectory, dct=None):
    """return metadata dictionary with defaults + custom dct"""

    attrs = {
        "System Name": get_sysname(trajectory.ref_atoms),
        "natoms": len(trajectory.ref_atoms),
        "time unit": "fs",
        "timestep": trajectory.timestep,
        "nsteps": len(trajectory) - 1,
        "volume": trajectory.volume,
        "symbols": trajectory.symbols,
        "masses": trajectory.masses,
        "flattened reference positions": trajectory.ref_positions.flatten(),
        "flattened lattice": trajectory.ref_atoms.cell.array.flatten(),
    }

    primitive = trajectory.primitive
    if primitive:
        prim_attrs = {
            "flattened primitive positions": primitive.positions.flatten(),
            "flattened primitive lattice": primitive.cell.array.flatten(),
        }
        attrs.update(prim_attrs)

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
        dims=vec_dims,
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
        dims=vec_dims,
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
        clean_pressure(trajectory.pressure / unit),
        dims=[time_dims],
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
        "forces": (vec_dims, trajectory.forces),
        "kinetic_energy": (time_dims, trajectory.kinetic_energy),
        "potential_energy": (time_dims, trajectory.potential_energy),
        "stress": (stress_dims, trajectory.stress),
        "pressure": pressure,
        "temperature": (time_dims, trajectory.temperatures),
    }
    coords = _time_coords(trajectory)
    attrs = _metadata(trajectory)

    return xr.Dataset(dataset, coords=coords, attrs=attrs)


def get_heat_flux_data(trajectory, return_avg=False, only_flux=False):
    """compute heat fluxes from TRAJECTORY and return as xarray

    Args:
        trajectory: list of atoms objects WITH ATOMIC STRESS computed
        return_avg (bool): return flux per atom
        only_flux (bool): only return heat flux and attrs

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
        "heat_flux": (vec_dims, trajectory.heat_flux),
        "pressure": data.pressure,
        "temperature": data.temperature,
    }

    if not only_flux:
        dataset.update(
            {
                "avg_heat_flux": (vec_dims, trajectory.avg_heat_flux),
                "positions": data.positions,
                "velocities": data.velocities,
                "forces": data.forces,
            }
        )
    coords = _time_coords(trajectory)
    attrs = _metadata(trajectory)

    return xr.Dataset(dataset, coords=coords, attrs=attrs)