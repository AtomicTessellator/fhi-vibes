#!/usr/bin/env python
"""Check the completion of an aims calculation"""
from ase.io import read

import balsam.launcher.dag as dag

from hilde.balsam.data_encoder import decode, encode
from hilde.balsam.post.aims_calc_process import postprocess_aims
from hilde.balsam.workflow.job_generator import generate_phonopy_jobs
from hilde.helpers.attribute_dict import AttributeDict
from hilde.phonon_db.ase_converters import dict2atoms, atoms2dict
from hilde.phonopy.context import PhonopyContext
from hilde.settings import Settings

completed, data = postprocess_aims(decode(dag.current_job.data))

if not completed:
    exit(0)

settings = Settings(read_config=False)
for key, val in data["ph_data"]["ctx"].items():
    if isinstance(val, dict):
        settings[key] = AttributeDict()
        for kk, vv in val.items():
            settings[key][kk] = vv
    else:
        settings[key] = val

settings["geometry"] = AttributeDict()
settings["geometry"]["file"] = "geometry.in"

atoms = dict2atoms(data["atoms"])
atoms.set_calculator(None)
atoms.write("geometry.in", scaled=True, geo_constrain=True)
atoms = read("geometry.in")

settings["control"].update(data["calculator"]["calculator_parameters"])
settings.pop("restart", None)
settings.pop("relaxation", None)
settings.general["relax_structure"] = False
settings.atoms = atoms

ctx = PhonopyContext(settings=settings, read_config=False)

supercell_calcs = generate_phonopy_jobs(ctx.settings, data)

for job in supercell_calcs:
    dag.add_dependency(dag.current_job, job)
    job_data = decode(job.data)
    job_data["primitive"] = atoms2dict(atoms)
    job.data = encode(job_data)

    job.save()
