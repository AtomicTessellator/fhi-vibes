#!/usr/bin/env python
"""Generate the aims.in file for hilde run aims"""

from balsam.launcher.dag import current_job

from hilde.aims.context import AimsContext
from hilde.aims.setup import setup_aims

from hilde.balsam.data_encoder import decode

from hilde.helpers.attribute_dict import AttributeDict
from hilde.phonon_db.ase_converters import dict2atoms
from hilde.settings import TaskSettings

data = decode(current_job.data)

atoms_dict = data["atoms"].copy()
for key, val in data["calculator"].items():
    atoms_dict[key] = val

settings = TaskSettings(read_config=False)
for key, val in data["ctx"].items():
    if isinstance(val, dict):
        settings[key] = AttributeDict()
        for kk, vv in val.items():
            settings[key][kk] = vv
    else:
        settings[key] = val

settings["control"].update(data["calculator"]["calculator_parameters"])
settings["geometry"] = AttributeDict({"file": "geometry.in"})
settings.write("aims.in")

atoms = dict2atoms(atoms_dict)
atoms.set_calculator(None)
atoms.write("geometry.in")

ctx = AimsContext(settings=settings)

calc = setup_aims(ctx)

calc.write_input(atoms)
