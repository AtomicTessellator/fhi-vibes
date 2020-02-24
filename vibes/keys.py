"""Keys and names for arrays. Composite keys are concatenated by `_`"""


def _join(*keys):
    """join two or more keys"""
    keys = (k for k in keys if k)
    return "_".join(keys)


fc = "force_constants"
fc_remapped = "force_constants_remapped"
fc_flattened = "force_constants_flattened"

system_name = "system_name"
time_unit = "time_unit"
timestep = "timestep"

reference_atoms = "atoms_reference"
reference_primitive = "atoms_primitive"
reference_supercell = "atoms_supercell"
reference_positions = "positions_reference"
reference_lattice = "lattice_reference"
metadata = "raw_metadata"

volume = "volume"
positions = "positions"
velocities = "velocities"

forces = "forces"
forces_harmonic = "forces_harmonic"
energy_kinetic = "energy_kinetic"
energy_potential = "energy_potential"
energy_potential_harmonic = "energy_potential_harmonic"
pressure = "pressure"
temperature = "temperature"

heat_flux = "heat_flux"
heat_fluxes = "heat_fluxes"
heat_flux_aux = "heat_flux_aux"
heat_fluxes_aux = "heat_fluxes_aux"

gk_prefactor = "gk_prefactor"
hfacf_scalar = "hfacf_scalar"
kappa_cumulative_scalar = "kappa_cumulative_scalar"
heat_flux_power_spectrum = "heat_flux_power_spectrum "

sigma = "sigma"
sigma_per_sample = "sigma_per_sample"

# time
time = "time"
omega = "omega"
cumtrapz = "cumtrapz"
autocorrelation = "autocorrelation"
fourier_transform = "fourier_transform"
avalanche_data = "avalanche_function"
time_avalanche = "avalanche_time"
avalanche_function = "avalanche_function"
avalanche_index = "avalanche_index"

# file management
default_backup_folder = "backups"
cache = "cache.vibes"

# hash
name = "name"
hash = "hash"
trajectory = "trajectory"


# composite keys
hfacf = _join(heat_flux, autocorrelation)
kappa_cumulative = _join(heat_flux, autocorrelation, cumtrapz)
