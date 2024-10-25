"""vibes Phono3py defaults"""

import collections

from vibes.helpers.dict import AttributeDict as adict

name = "phono3py"

mandatory = {
    "mandatory_keys": ["machine", "control", "geometry"],
    "mandatory_obj_keys": ["supercell_matrix"],
}

_keys = [
    "supercell_matrix",
    "displacement",
    "cutoff_pair_distance",
    "symprec",
    "is_diagonal",
    "is_plusminus",
    "log_level",
    "workdir",
]
keys = collections.namedtuple(f"{name}_keywords", _keys)(*_keys)

kwargs = adict(
    {
        keys.supercell_matrix: [1, 1, 1],
        keys.displacement: 0.03,
        keys.cutoff_pair_distance: 100.0,
        keys.symprec: 1e-5,
        keys.is_diagonal: True,
        keys.is_plusminus: "auto",
        keys.log_level: 2,
        keys.workdir: name,
    }
)

fc2_file = "fc2.hdf5"
fc3_file = "fc3.hdf5"
disp_fc3_yaml_file = "disp_fc3.yaml"
phono3py_params_yaml_file = "phono3py_params.yaml"

settings_dict = {name: kwargs}
