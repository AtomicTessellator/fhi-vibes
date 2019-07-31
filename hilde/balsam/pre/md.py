#!/usr/bin/env python3
'''Generate the aims.in file for hilde run aims'''

import numpy as np
from balsam.launcher.dag import current_job

from hilde.helpers.attribute_dict import AttributeDict
from hilde.phonon_db.ase_converters import dict2atoms
from hilde.settings import TaskSettings

data = current_job.data

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

set.write("md.in")
