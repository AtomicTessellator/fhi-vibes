""" hilde defaults for aims"""

from hilde.helpers.attribute_dict import AttributeDict as adict

name = "aims"

obj_key = "control"
mandatory_base = ["machine", "control", "geometry"]
mandatory_task = ["xc"]

# either of
mandatory_basisset = ("basisset", "basissets")


defaults = adict(
    {
        "sc_accuracy_rho": 1e-6,
        "relativistic": "atomic_zora scalar",
        "output_level": "MD_light",
    }
)
