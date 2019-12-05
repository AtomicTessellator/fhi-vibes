""" vibes defaults for md"""

from vibes.helpers.attribute_dict import AttributeDict as adict

name = "md"

mandatory_base = ["machine", "geometry", name]
mandatory_task = ["driver", "timestep", "maxsteps"]

defaults = adict({"driver": "VelocityVerlet", "logfile": "md.log"})
