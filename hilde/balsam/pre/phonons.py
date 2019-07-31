#!/usr/bin/env python
'''Generate the aims.in file for hilde run aims'''

import numpy as np
from balsam.launcher.dag import current_job

from hilde.balsam.data_encoder import decode

from hilde.helpers.attribute_dict import AttributeDict
from hilde.phonon_db.ase_converters import dict2atoms
from hilde.settings import TaskSettings

data = decode(current_job.data)

atoms_dict = data["atoms"].copy()
for key, val in data["calculator"].items():
    atoms_dict[key] = val

set = TaskSettings(read_config=False)
for key, val in data["ctx"].items():
    if isinstance(val, dict):
        set[key] = AttributeDict()
        for kk, vv in val.items():
            set[key][kk] = vv
    else:
        set[key] = val

set["geometry"] = AttributeDict()
set["geometry"]["file"] = "geometry.in"
dict2atoms(data["atoms"]).write("geometry.in", scaled=True, use_sym=True)

set["control"].update(data["calculator"]["calculator_parameters"])

if "restart" in set:
    del(set["restart"])

set.write("phonopy.in")
