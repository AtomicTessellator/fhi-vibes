""" tools for storing MD trajectories

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

# from vibes import __version__ as version
# from vibes import son
# from vibes.helpers.converters import results2dict
# from vibes.helpers.converters import dict2json as dumper
from vibes.helpers import talk as _talk, Timer

_prefix = "trajectory"
_fc_key = "force_constants"

key_reference_atoms = "reference atoms"
key_reference_primitive = "reference primitive atoms"
key_reference_positions = "reference positions"
key_metadata = "raw_metadata"

key_positions = "positions"
key_velocities = "velocities"

key_forces = "forces"
key_forces_harmonic = "forces_harmonic"
key_energy_kinetic = "kinetic_energy"
key_energy_potential = "potential_energy"
key_energy_potential_harmonic = "potential_energy_harmonic"
key_pressure = "pressure"
key_temperature = "temperature"

time_dims = "time"
atoms_dims = [time_dims, "I"]
vec_dims = [time_dims, "a"]
atoms_vec_dims = [time_dims, "I", "a"]
stress_dims = [time_dims, "a", "b"]
stresses_dims = [time_dims, "I", "a", "b"]
kappa_dims = [time_dims, "I", "J", "a", "b"]


Timer.prefix = _prefix


def talk(msg):
    """wrapper for `utils.talk` with prefix"""
    return _talk(msg, prefix=_prefix)


def get_hashes_from_trajectory_file(trajectory, verbose=False):
    """return all hashes from trajectory"""

    try:
        traj = reader(trajectory, verbose=verbose)
    except (FileNotFoundError, KeyError):
        return []

    return traj.get_hashes()


from vibes.helpers.hash import hash_atoms
from vibes.trajectory.io import reader, metadata2file, step2file, results2dict
from vibes.trajectory.trajectory import Trajectory
