#!/usr/bin/env python
"""Check the completion of an aims calculation"""


from glob import glob

from pathlib import Path

from shutil import copyfile

import balsam.launcher.dag as dag

import numpy as np

from hilde.aims.context import AimsContext
from hilde.aims.setup import setup_aims
from hilde.balsam.data_encoder import decode
from hilde.balsam.post.aims_calc_process import postprocess_aims
from hilde.balsam.post.phonopy_postprocess_functions import (
    collect_calcs_to_trajectory,
    setup_gruneisen,
    setup_statistical_sampling,
)
from hilde.balsam.workflow.job_generator import (
    generate_phonopy_jobs,
    get_phonon_setup_data,
)
from hilde.fireworks.tasks.postprocess.calculate import get_calc_times
from hilde.fireworks.tasks.postprocess.phonons import get_converge_phonon_update
from hilde.fireworks.workflows.firework_generator import time2str
from hilde.helpers.attribute_dict import AttributeDict
from hilde.helpers.paths import cwd
from hilde.phonon_db.ase_converters import dict2atoms
from hilde.phonopy.context import PhonopyContext
from hilde.phonopy.postprocess import postprocess
from hilde.settings import Settings

data = decode(dag.current_job.data)

# Postprocess the aims calculation
completed, _ = postprocess_aims(decode(dag.current_job.data))

if not completed:
    exit(0)

# Move
Path(data["analysis_wd"]).mkdir(exist_ok=True, parents=True)
copyfile("aims.out", f"{data['analysis_wd']}/aims.out")

base_dir = data["ph_data"]["analysis_wd"]
completed_calcs = glob(f"{base_dir}/*/aims.out")
if len(completed_calcs) < data["ph_data"]["number_of_sc"]:
    exit(0)

calc_dirs = glob(f"{base_dir}/*/")

with cwd(base_dir):
    collect_calcs_to_trajectory(base_dir, data["ph_data"]["metadata"])
    ctx = PhonopyContext(
        Settings(settings_file="phonopy.in", read_config=False), read_config=False
    )
    calc = setup_aims(AimsContext(Settings(settings_file="phonopy.in")))
    ctx.settings.atoms.set_calculator(calc)

    atoms = ctx.settings.atoms
    atoms.set_calculator(calc)
    atoms.info.update(data["atoms"]["info"])

    workdir = data["ph_data"]["analysis_wd"]
    trajectory = "trajectory.son"
    ph_times = get_calc_times(calc_dirs=calc_dirs)

    phonon = postprocess(f"./{trajectory}")

    conv, update_job = get_converge_phonon_update(
        workdir,
        trajectory,
        ph_times,
        phonon,
        data["ph_data"].get("conv_crit", None),
        data["ph_data"].get("prev_dos_fp", None),
        sc_matrix_original=data["ph_data"].get("sc_matrix_original", None),
    )

    update_job["ph_calculator"]["calculator_parameters"].pop("k_grid", None)
    update_job["ph_calculator"]["calculator_parameters"].pop("use_pimd_wrapper", None)

    update_data = {
        "calculator": update_job["ph_calculator"],
        "k_pt_density": update_job["k_pt_density"],
    }

    if conv or not ctx.settings.phonopy.converge_phonons:
        if ctx.settings.phonopy.get("converge_phonons", False):
            sc_matrix_max = np.max(ctx.settings.phonopy.supercell_matrix.flatten())
            sc_matrix_orig_max = np.max(data["ph_data"]["sc_matrix_original"].flatten())
            n = int(round(sc_matrix_max / sc_matrix_orig_max))
            n -= 1
            sc_mat = n * np.array(data["ph_data"]["sc_matrix_original"]).reshape((3, 3))
        else:
            sc_mat = ctx.settings.phonopy.supercell_matrix

        if ctx.settings.phonopy.get_gruniesen:
            with cwd(base_dir):
                setup_gruneisen(0.99, sc_mat)
                setup_gruneisen(1.01, sc_mat)

        if "statistical_sampling" in ctx.settings:
            with cwd(base_dir):
                setup_statistical_sampling(sc_mat, phonon)
        exit(0)

update_job["expected_walltime"] /= data["ph_data"]["number_of_sc"]

ctx.settings.atoms = dict2atoms(data["primitive"])
ctx.settings.control.update(update_data["calculator"]["calculator_parameters"])
ctx.settings["control_kpt"] = AttributeDict(dict(density=update_data["k_pt_density"]))
ctx.ref_atoms = dict2atoms(data["primitive"])

qadapter = dict(walltime=time2str(update_job["expected_walltime"]))
ctx.settings.phonopy["supercell_matrix"] = update_job["supercell_matrix"]

ph_data = get_phonon_setup_data(ctx.settings, qadapter, update_job["prev_dos_fp"])
data.update(update_data)
data["ph_data"].update(ph_data)
jobs = generate_phonopy_jobs(ctx.settings, data)

for job in jobs:
    dag.add_dependency(dag.current_job, job)
    job.wall_time_minutes = float(dag.current_job.num_nodes) / float(job.num_nodes)
    job.save()
