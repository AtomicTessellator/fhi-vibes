""" hilde defaults for md"""

from hilde.helpers.attribute_dict import AttributeDict as adict

name = "md"

mandatory_base = ["machine", "control", "geometry", name]
mandatory_task = ["driver", "timestep", "maxsteps"]

defaults = adict({"driver": "VelocityVerlet", "logfile": "md.log"})
