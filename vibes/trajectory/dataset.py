"""compute and analyze heat fluxes"""
import json
import numpy as np
import xarray as xr

from ase import units
from vibes.helpers import warn
from vibes.helpers.converters import atoms2json
from vibes.structure.misc import get_sysname
from .utils import clean_pressure
from vibes.green_kubo.heat_flux import (
    key_heat_flux,
    key_heat_flux_aux,
    key_heat_fluxes,
    key_heat_fluxes_aux,
)
from . import (
    Timer,
    key_reference_atoms,
    key_reference_positions,
    key_reference_primitive,
    key_metadata,
    key_forces,
    key_forces_harmonic,
    key_energy_kinetic,
    key_energy_potential,
    key_energy_potential_harmonic,
    time_dims,
    vec_dims,
    atoms_vec_dims,
    stress_dims,
    kappa_dims,
)


def _time_coords(trajectory):
    """return time as coords dict"""
    coords = {time_dims: trajectory.times}
    return coords


def _attrs(trajectory, dct=None, metadata=False):
    """return metadata dictionary with defaults + custom dct"""

    attrs = {
        "System Name": get_sysname(trajectory.ref_atoms),
        "natoms": len(trajectory.ref_atoms),
        "time unit": "fs",
        "timestep": trajectory.timestep,
        "nsteps": len(trajectory) - 1,
        "symbols": trajectory.symbols,
        "masses": trajectory.masses,
        key_reference_atoms: atoms2json(trajectory.reference_atoms, reduce=False),
        "average atoms": atoms2json(trajectory.average_atoms, reduce=False),
        key_reference_positions: trajectory.ref_positions.flatten(),
    }

    # handle non-periodic systems
    try:
        attrs.update({"volume": trajectory.volume})
    except ValueError:
        pass

    if trajectory.primitive:
        rep = atoms2json(trajectory.primitive, reduce=False)
        prim_attrs = {key_reference_primitive: rep}
        attrs.update(prim_attrs)

    if trajectory.force_constants_remapped is not None:
        fc = trajectory.force_constants_remapped.flatten()
        fc_attrs = {"flattened force_constants": fc}
        attrs.update(fc_attrs)

    if dct and isinstance(dct, dict):
        attrs.update(dct)

    if metadata:
        raw_metadata = json.dumps(trajectory.metadata)
        attrs.update({key_metadata: raw_metadata})

    return attrs


def get_positions_dataarray(trajectory, verbose=True):
    """extract positions from TRAJECTORY  and return as xarray.DataArray

    Args:
        trajectory (Trajectory): list of atoms objects
    Returns:
        positions (xarray.DataArray [N_t, N_a, 3])
    """
    timer = Timer("Get positions from trajectory", verbose=verbose)

    df = xr.DataArray(
        trajectory.positions,
        dims=atoms_vec_dims,
        coords=_time_coords(trajectory),
        name="positions",
        attrs=_attrs(trajectory),
    )

    timer()

    return df


def get_velocities_dataarray(trajectory, verbose=True):
    """extract velocties from TRAJECTORY  and return as xarray.DataArray

    Args:
        trajectory (Trajectory): list of atoms objects
    Returns:
        velocities (xarray.DataArray [N_t, N_a, 3])
    """
    timer = Timer("Get velocities from trajectory", verbose=verbose)

    df = xr.DataArray(
        trajectory.velocities,
        dims=atoms_vec_dims,
        coords=_time_coords(trajectory),
        name="velocities",
        attrs=_attrs(trajectory),
    )

    timer()

    return df


def get_pressure_dataarray(trajectory, GPa=False, verbose=True):
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
        attrs=_attrs(trajectory, dct=extra_metadata),
    )

    timer()

    return df


def get_trajectory_dataset(trajectory, metadata=False):
    """Return trajectory data as xarray.Dataset

    Args:
        trajectory: list of atoms objects WITH ATOMIC STRESS computed
        metadata (bool): include `raw_metadata` in `attrs`
    Returns:
        xarray.Dataset:
            positions, velocities, forces, stress, pressure, temperature
    """

    # add velocities and pressure
    positions = get_positions_dataarray(trajectory)
    velocities = get_velocities_dataarray(trajectory)
    pressure = get_pressure_dataarray(trajectory)

    dataset = {
        "positions": positions,
        "displacements": (atoms_vec_dims, trajectory.displacements),
        "velocities": velocities,
        "momenta": (atoms_vec_dims, trajectory.momenta),
        key_forces: (atoms_vec_dims, trajectory.forces),
        key_energy_kinetic: (time_dims, trajectory.kinetic_energy),
        key_energy_potential: (time_dims, trajectory.potential_energy),
        "stress": (stress_dims, trajectory.stress),
        "pressure": pressure,
        "temperature": (time_dims, trajectory.temperatures),
    }

    coords = _time_coords(trajectory)
    attrs = _attrs(trajectory, metadata=metadata)

    if trajectory.forces_harmonic is not None:
        epot_ha = trajectory.potential_energy_harmonic
        update_dict = {
            key_forces_harmonic: (atoms_vec_dims, trajectory.forces_harmonic),
            key_energy_potential_harmonic: (time_dims, epot_ha),
            "sigma_per_sample": (time_dims, trajectory.sigma_per_sample),
        }
        dataset.update(update_dict)
        attrs.update({"sigma": trajectory.sigma})

    return xr.Dataset(dataset, coords=coords, attrs=attrs)


def get_heat_flux_dataset(trajectory, only_flux=False):
    """compute heat fluxes from TRAJECTORY and return as xarray

    Args:
        trajectory: list of atoms objects WITH ATOMIC STRESS computed
        only_flux (bool): only return heat flux and attrs

    Returns:
        xarray.Dataset:
            heat_flux
            avg_heat_flux
    """
    # add velocities and pressure
    data = get_trajectory_dataset(trajectory)

    flux = [a.calc.results[key_heat_flux] for a in trajectory]

    dataset = {
        key_heat_flux: (vec_dims, np.array(flux)),
        "pressure": data.pressure,
        "temperature": data.temperature,
    }

    if not only_flux:
        fluxes = [a.calc.results[key_heat_fluxes] for a in trajectory]
        flux_aux = [a.calc.results[key_heat_flux_aux] for a in trajectory]
        fluxes_aux = [a.calc.results[key_heat_fluxes_aux] for a in trajectory]

        dataset.update(
            {
                key_heat_fluxes: (atoms_vec_dims, np.array(fluxes)),
                key_heat_flux_aux: (vec_dims, np.array(flux_aux)),
                key_heat_fluxes_aux: (atoms_vec_dims, np.array(fluxes_aux)),
                "positions": data.positions,
                "velocities": data.velocities,
                key_forces: data.forces,
            }
        )
    coords = _time_coords(trajectory)
    attrs = _attrs(trajectory)

    return xr.Dataset(dataset, coords=coords, attrs=attrs)
