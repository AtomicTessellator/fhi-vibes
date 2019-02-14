""" hilde defaults for phonopy """

from hilde.helpers.attribute_dict import AttributeDict as adict

displacement_id_str = "displacement_id"

defaults = adict(
    {
        # for phono3py compatibility
        "displacement": 0.03,
        "symprec": 1e-5,
        "trigonal": False,
        "is_diagonal": False,
        "q_mesh": [11, 11, 11],
    }
)