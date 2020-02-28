""" vibes Phono3py defaults """

from vibes.phonopy import defaults as defaults_fc2

name = "phono3py"

mandatory = {
    "mandatory_keys": ["machine", "control", "geometry"],
    "mandatory_obj_keys": ["supercell_matrix"],
}

defaults = defaults_fc2.copy()
defaults.update(
    {
        "displacement": 0.03,
        "is_diagonal": True,
        "cutoff_pair_distance": 100.0,
        "qmesh": [21, 21, 21],
        "log_level": 2,
    }
)
