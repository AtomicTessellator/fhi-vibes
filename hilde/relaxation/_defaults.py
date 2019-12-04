""" hilde defaults for md"""

from hilde.helpers.attribute_dict import AttributeDict as adict

name = "relaxation"

mandatory_base = ["machine", "control", "geometry", name]
mandatory_task = ["driver", "fmax"]

defaults = adict(
    {
        "driver": "BFGS",
        "logfile": "relaxation.log",
        "unit_cell": False,
        "fmax": 0.001,
        # "alpha": 25,
        "maxstep": 0.04,
        "decimals": 10,
    }
)
