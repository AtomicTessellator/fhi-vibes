""" vibes defaults for aims"""

from vibes.helpers import talk as _talk
from vibes.helpers.dict import AttributeDict as adict

name = "aims"

obj_key = "control"
mandatory_base = ["machine", "control", "geometry"]
mandatory_task = ["xc"]

basisset_key = "basissets"
basisset_choices = ("light", "intermediate", "tight", "really_tight")
basisset_default = "light"


kwargs = adict(
    {
        "sc_accuracy_rho": 1e-6,
        "relativistic": "atomic_zora scalar",
        "output_level": "MD_light",
    }
)


def talk(msg, verbose=True):
    """wrapper for helpers.talk with 'aims' prefix"""
    return _talk(msg, prefix=name, verbose=verbose)
