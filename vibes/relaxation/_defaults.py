""" vibes defaults for md"""

from vibes.helpers.attribute_dict import AttributeDict as adict

name = "relaxation"

mandatory_base = ["machine", "geometry", name]
mandatory_task = ["driver", "fmax"]

kwargs= adict(
    {
        "driver": "BFGS",
        "logfile": "relaxation.log",
        "unit_cell": False,
        "fmax": 0.001,
        # "alpha": 25,
        "maxstep": 0.04,
        "hydrostatic_strain": False,
        "constant_volume": False,
        "scalar_pressure": 0.0,
        "decimals": 10,
    }
)
