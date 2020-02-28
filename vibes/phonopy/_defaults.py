""" vibes defaults for phonopy """

from vibes.helpers.attribute_dict import AttributeDict as adict

displacement_id_str = "displacement_id"
name = "phonopy"

mandatory = {
    "mandatory_keys": ["machine", "control", "geometry"],
    "mandatory_obj_keys": ["supercell_matrix"],
}

defaults = adict(
    {
        # for phono3py compatibility
        "displacement": 0.01,
        "symprec": 1e-5,
        "is_diagonal": False,
        "is_trigonal": False,
        "is_plusminus": "auto",
        "q_mesh": [45, 45, 45],
    }
)
