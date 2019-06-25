""" hilde defaults for phonopy """

from hilde.helpers.attribute_dict import AttributeDict as adict

displacement_id_str = "displacement_id"
name = "phonopy"

mandatory_base = ["machine", "control", "geometry", name]
mandatory_task = ["supercell_matrix"]

defaults = adict(
    {
        # for phono3py compatibility
        "displacement": 0.01,
        "symprec": 1e-5,
        "is_diagonal": False,
        "q_mesh": [26, 26, 26],
    }
)
