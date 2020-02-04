"""compute and analyze heat fluxes"""
import json

import numpy as np
import xarray as xr
from ase import units

from vibes import dimensions as dims
from vibes import keys
from vibes.helpers import warn
from vibes.helpers.converters import atoms2json, dict2json
from vibes.structure.misc import get_sysname

from .utils import Timer, clean_pressure


def _time_coords(trajectory):
    """return time as coords dict"""
    coords = {dims.time: trajectory.times}
    return coords


def _attrs(trajectory, dct=None, metadata=False):
    """return metadata dictionary with defaults + custom dct"""

    attrs = {
        keys.system_name: get_sysname(trajectory.ref_atoms),
        "natoms": len(trajectory.ref_atoms),
        keys.time_unit: "fs",
        keys.timestep: trajectory.timestep,
        "nsteps": len(trajectory) - 1,
        "symbols": trajectory.symbols,
        "masses": trajectory.masses,
        keys.reference_atoms: atoms2json(trajectory.reference_atoms, reduce=False),
    }

    if trajectory.primitive:
        rep = atoms2json(trajectory.primitive, reduce=False)
        prim_attrs = {keys.reference_primitive: rep}
        attrs.update(prim_attrs)

    if trajectory.supercell:
        rep = atoms2json(trajectory.supercell, reduce=False)
        prim_attrs = {keys.reference_supercell: rep}
        attrs.update(prim_attrs)

    # handle non-periodic systems
    try:
        attrs.update({"volume": trajectory.volume})
    except ValueError:
        pass

    if dct and isinstance(dct, dict):
        attrs.update(dct)

    if metadata:
        raw_metadata = dict2json(trajectory.metadata)
        attrs.update({keys.metadata: raw_metadata})

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
        dims=dims.atoms_vec,
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
        dims=dims.atoms_vec,
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
        dims=[dims.time],
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

    # reference positions
    positions_reference = (dims.positions, trajectory.reference_atoms.positions)
    lat = np.asarray(trajectory.reference_atoms.cell)
    lattice_reference = (dims.lattice, lat)

    dataset = {
        "positions": positions,
        "displacements": (dims.atoms_vec, trajectory.displacements),
        "velocities": velocities,
        "momenta": (dims.atoms_vec, trajectory.momenta),
        keys.forces: (dims.atoms_vec, trajectory.forces),
        keys.energy_kinetic: (dims.time, trajectory.kinetic_energy),
        keys.energy_potential: (dims.time, trajectory.potential_energy),
        "stress": (dims.stress, trajectory.stress),
        "pressure": pressure,
        "temperature": (dims.time, trajectory.temperatures),
        keys.reference_positions: positions_reference,
        keys.reference_lattice: lattice_reference,
    }

    # heat_flux
    flux = trajectory.get_heat_flux()
    if flux is not None:
        dataset.update({keys.heat_flux: (dims.vec, flux)})

    # heat_flux_aux
    flux = trajectory.get_heat_flux(aux=True)
    if flux is not None:
        dataset.update({keys.heat_flux_aux: (dims.vec, flux)})

    # heat_fluxes
    flux = trajectory.get_heat_fluxes()
    if flux is not None:
        dataset.update({keys.heat_fluxes: (dims.atoms_vec, flux)})

    # heat_fluxes_aux
    flux = trajectory.get_heat_fluxes(aux=True)
    if flux is not None:
        dataset.update({keys.heat_fluxes_aux: (dims.atoms_vec, flux)})

    coords = _time_coords(trajectory)
    attrs = _attrs(trajectory, metadata=metadata)

    if trajectory.force_constants_remapped is not None:
        fc = trajectory.force_constants_remapped
        dataset.update({keys.fc_remapped: (dims.fc_remapped, fc)})

    if trajectory.forces_harmonic is not None:
        epot_ha = trajectory.potential_energy_harmonic
        update_dict = {
            keys.forces_harmonic: (dims.atoms_vec, trajectory.forces_harmonic),
            keys.energy_potential_harmonic: (dims.time, epot_ha),
            keys.sigma_per_sample: (dims.time, trajectory.sigma_per_sample),
        }
        dataset.update(update_dict)
        attrs.update({"sigma": trajectory.sigma})

    return xr.Dataset(dataset, coords=coords, attrs=attrs)


def get_heat_flux_dataset(trajectory, only_flux=False, metadata=False):
    """compute heat fluxes from TRAJECTORY and return as xarray

    Args:
        trajectory: list of atoms objects WITH ATOMIC STRESS computed
        only_flux (bool): only return heat flux and attrs
        metadata (bool): include `raw_metadata` in `attrs`

    Returns:
        xarray.Dataset:
            heat_flux
            avg_heat_flux
    """
    # add velocities and pressure
    data = get_trajectory_dataset(trajectory, metadata=metadata)

    flux = [a.calc.results[keys.heat_flux] for a in trajectory]

    dataset = {
        keys.heat_flux: (dims.vec, np.array(flux)),
        "pressure": data.pressure,
        "temperature": data.temperature,
    }

    if not only_flux:
        fluxes = [a.calc.results[keys.heat_fluxes] for a in trajectory]
        flux_aux = [a.calc.results[keys.heat_flux_aux] for a in trajectory]
        fluxes_aux = [a.calc.results[keys.heat_fluxes_aux] for a in trajectory]

        dataset.update(
            {
                keys.heat_fluxes: (dims.atoms_vec, np.array(fluxes)),
                keys.heat_flux_aux: (dims.vec, np.array(flux_aux)),
                keys.heat_fluxes_aux: (dims.atoms_vec, np.array(fluxes_aux)),
                "positions": data.positions,
                "velocities": data.velocities,
                keys.forces: data.forces,
            }
        )
    coords = _time_coords(trajectory)
    attrs = _attrs(trajectory)

    return xr.Dataset(dataset, coords=coords, attrs=attrs)
