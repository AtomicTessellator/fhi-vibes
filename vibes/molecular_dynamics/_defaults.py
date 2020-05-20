""" vibes defaults for md"""
import collections

from vibes.helpers.dict import AttributeDict as adict

name = "md"

_keys = [
    "driver",
    "timestep",
    "temperature",
    "maxsteps",
    "compute_stresses",
    "logfile",
    "friction",
    "workdir",
]
keys = collections.namedtuple("md_keywords", _keys)(*_keys)

kwargs = adict(
    {
        keys.driver: "VelocityVerlet",
        keys.timestep: 1,
        keys.maxsteps: 1000,
        keys.compute_stresses: False,
        keys.workdir: name,
        # kwargs go to Dynamics, e.g., Langeving(..., **kwargs)
        "kwargs": {keys.temperature: None, keys.friction: None, keys.logfile: "md.log"},
    }
)

settings_dict = {name: kwargs}

calculation_timeout = 30 * 60  # 30 minutes
