"""compute and analyze heat fluxes"""
import numpy as np
import xarray as xr
from vibes import dimensions as dims
from vibes import keys
from vibes.helpers import warn
from vibes.helpers.converters import atoms2json, dict2json
from vibes.structure.misc import get_sysname

from .utils import Timer


def _time_coords(trajectory):
    """return time as coords dict"""
    coords = {dims.time: trajectory.times}
    return coords


def _attrs(trajectory, dct=None, metadata=False):
    """return metadata dictionary with defaults + custom dct"""

    attrs = {
        keys.name: keys.trajectory,
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
        attrs.update({keys.volume: trajectory.volume})
    except ValueError:
        pass

    if dct and isinstance(dct, dict):
        attrs.update(dct)

    if metadata:
        raw_metadata = dict2json(trajectory.metadata)
        attrs.update({keys.metadata: raw_metadata})

    # is slow:
    # attrs.update({keys.hash: trajectory.hash})  # add hash
    if trajectory.hash_raw:
        attrs.update({keys.hash_raw: trajectory.hash_raw})  # add raw hash

    return attrs


def get_pressure_dataset(trajectory, attrs=None, verbose=True):
    """extract pressure from TRAJECTORY  and return as xarray.DataArray

    Args:
        trajectory (Trajectory): list of atoms objects
        attrs: attach attrs (if already obtained)
    Returns:
        pressure (xarray.DataArray [N_t]) in eV/AA**3
    """
    timer = Timer("Get pressure from trajectory", verbose=verbose)

    data = {
        keys.pressure: (dims.time, trajectory.pressure),
        keys.pressure_kinetic: (dims.time, trajectory.pressure_kinetic),
        keys.pressure_potential: (dims.time, trajectory.pressure_potential),
    }

    if attrs is None:
        attrs = _attrs(trajectory)

    df = xr.Dataset(
        data,
        coords=_time_coords(trajectory),
        attrs=attrs,
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

    # reference positions
    positions_reference = (dims.positions, trajectory.reference_atoms.positions)
    lat = np.asarray(trajectory.reference_atoms.cell)
    lattice_reference = (dims.lattice, lat)

    dataset = {
        keys.reference_positions: positions_reference,
        keys.reference_lattice: lattice_reference,
        keys.positions: (dims.time_atom_vec, trajectory.positions),
        keys.displacements: (dims.time_atom_vec, trajectory.displacements),
        keys.velocities: (dims.time_atom_vec, trajectory.velocities),
        keys.momenta: (dims.time_atom_vec, trajectory.momenta),
        keys.forces: (dims.time_atom_vec, trajectory.forces),
        keys.energy_kinetic: (dims.time, trajectory.kinetic_energy),
        keys.energy_potential: (dims.time, trajectory.potential_energy),
        keys.temperature: (dims.time, trajectory.temperatures),
        keys.cell: (dims.time_tensor, trajectory.cells),
        keys.volume: (dims.time, trajectory.volumes),
        keys.stress: (dims.time_tensor, trajectory.stress),
        keys.stress_potential: (dims.time_tensor, trajectory.stress_potential),
        keys.stress_kinetic: (dims.time_tensor, trajectory.stress_kinetic),
    }

    stresses_potential = trajectory.stresses_potential
    if stresses_potential is not None:
        value = (dims.time_atom_tensor, trajectory.stresses_potential)
        dataset.update({keys.stresses_potential: value})

    virials = trajectory.virials
    if virials is not None:
        value = (dims.time_atom_tensor, trajectory.virials)
        dataset.update({keys.virials: value})

    # heat_flux
    flux = trajectory.get_heat_flux()
    if flux is not None:
        dataset.update({keys.heat_flux: (dims.time_vec, flux)})

    # heat_flux_aux
    flux = trajectory.get_heat_flux(aux=True)
    if flux is not None:
        dataset.update({keys.heat_flux_aux: (dims.time_vec, flux)})

    # heat_flux_harmonic
    flux = trajectory.get_heat_flux_harmonic()
    if flux is not None:
        dataset.update({keys.heat_flux_harmonic: (dims.time_vec, flux)})
    flux = trajectory.get_heat_flux_0_harmonic()
    if flux is not None:
        dataset.update({keys.heat_flux_0_harmonic: (dims.time_vec, flux)})

    coords = _time_coords(trajectory)
    attrs = _attrs(trajectory, metadata=metadata)

    if trajectory.force_constants is not None:
        fc = trajectory.force_constants.array
        dataset.update({keys.fc: (dims.fc, fc)})
        rfc = trajectory.force_constants_remapped
        dataset.update({keys.fc_remapped: (dims.fc_remapped, rfc)})
        map_s2p = trajectory.force_constants.I2iL_map[:, 0]
        attrs.update({keys.map_supercell_to_primitive: map_s2p})

    if trajectory.forces_harmonic is not None:
        epot_ha = trajectory.potential_energy_harmonic
        update_dict = {
            keys.forces_harmonic: (dims.time_atom_vec, trajectory.forces_harmonic),
            keys.energy_potential_harmonic: (dims.time, epot_ha),
            keys.sigma_per_sample: (dims.time, trajectory.sigma_per_sample),
        }
        dataset.update(update_dict)
        attrs.update({keys.sigma: trajectory.sigma})

    ds = xr.Dataset(dataset, coords=coords, attrs=attrs)

    # add pressure
    ds.update(get_pressure_dataset(trajectory, attrs=attrs))

    # aims uuid
    aims_uuids = trajectory.aims_uuid
    if aims_uuids[0] is not None:
        ds.update({keys.aims_uuid: (dims.time, aims_uuids)})

    # check for duplicate times
    unique_times, index = np.unique(ds[keys.time], return_index=True)
    dn = len(ds[keys.time]) - len(unique_times)
    if dn > 0:
        warn(f"Drop {dn} duplicate timesteps. CHECK TRAJECTORY!", level=1)
        ds = ds.isel({keys.time: index})

    return ds
