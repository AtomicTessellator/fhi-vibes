#!/usr/bin/env python
"""Check the completion of an aims calculation"""

import balsam.launcher.dag as dag

import numpy as np

from hilde.aims.context import AimsContext
from hilde.aims.setup import setup_aims

from hilde.balsam.data_encoder import decode, encode
from hilde.balsam.workflow.job_generator import (
    generate_gruneisen_jobs,
    generate_phonopy_job,
)
from hilde.fireworks.tasks.postprocess.calculate import get_calc_times
from hilde.fireworks.tasks.postprocess.phonons import get_converge_phonon_update
from hilde.materials_fp import material_fingerprint
from hilde.phonopy.context import PhonopyContext
from hilde.phonopy.postprocess import postprocess
from hilde.phonopy.workflow import bootstrap
from hilde.settings import Settings
from hilde.trajectory import reader


def setup_gruneisen(vol_factor, sc_mat):
    ctx = PhonopyContext(Settings(settings_file="phonopy.in", read_config=False))
    calc = setup_aims(AimsContext(Settings(settings_file="phonopy.in")))
    ctx.settings.atoms.set_calculator(calc)

    ctx.settings.phonopy.supercell_matrix = sc_mat
    ctx.settings.phonopy["get_gruniesen"] = False
    ctx.settings.phonopy["converge_phonons"] = False

    jobs, job_dep = generate_gruneisen_jobs(ctx.settings, vol_factor)

    for job in jobs:
        job.save()

    for key, val in job_dep.items():
        dag.add_dependency(key, val)

    dag.add_dependency(dag.current_job, jobs[0])
    for child in dag.children:
        dag.add_dependency(jobs[-1], child)

    return jobs, job_dep


data = decode(dag.current_job.data)

ctx = PhonopyContext(Settings(settings_file="phonopy.in"))
calc = setup_aims(AimsContext(Settings(settings_file="phonopy.in")))
ctx.settings.atoms.set_calculator(calc)

atoms = ctx.settings.get_atoms()
atoms.set_calculator(calc)
atoms.info.update(data["atoms"]["info"])

workdir = "./phonopy/"
trajectory = "trajectory.son"
ph_times = get_calc_times(workdir)

try:
    ph = postprocess(f"{workdir}/{trajectory}")
except RuntimeError:
    ca = reader("trajectory.son")
    atoms_to_calculate = bootstrap(ctx)["atoms_to_calculate"]
    if len(ca) > 0:
        scale_wt = (len(atoms_to_calculate) - len(ca)) / len(ca)
    else:
        raise RuntimeError("Force calculation ran out of time for first calculation.")
    for child in dag.children:
        dag.kill(child, recursive=True)
    dag.spawn_child(
        clone=True, walltime_minutes=dag.current_job.walltime_minutes * scale_wt
    )

conv, update_job = get_converge_phonon_update(
    workdir,
    trajectory,
    ph_times,
    ph,
    data["conv_crit"],
    data["prev_dos_fp"],
    workdir,
    sc_matrix_original=data["sc_matrix_original"],
)
update_data = {
    "atoms": update_job["ph_primitive"],
    "calculator": update_job["ph_calculator"],
    "k_pt_density": update_job["k_pt_density"],
}
update_data["atoms"]["info"] = data["atoms"]["info"]
if not conv:
    update_data["prev_dos_fp"] = update_job["prev_dos_fp"]

for child in dag.children:
    ch_data = decode(child.data)
    ch_data.update(update_job)
    child.data = encode(ch_data.copy())
    child.save()

if conv or not ctx.settings.phonopy.converge_phonons:
    if ctx.settings.phonopy.converge_phonons:
        update_data["ph_time"] = update_job["ph_time"]
    if getattr(ctx.settings.phonopy, "get_gruniesen", False):
        if ctx.settings.phonopy.converge_phonons:
            n = ph.get_supercell_matrix()[0, 0] / data["sc_matrix_original"][0]
            n -= 1
            sc_mat = n * np.array(data["sc_matrix_original"]).reshape((3, 3))
        else:
            sc_mat = ph.get_supercell_matrix()
        jobs_mn, dep_mn = setup_gruneisen(0.99, sc_mat)
        jobs_pl, dep_pl = setup_gruneisen(1.01, sc_mat)
    exit(0)

ctx.settings.phonopy.supercell_matrix = list(update_job["supercell_matrix"].flatten())
natoms = len(update_job["ph_supercell"]["numbers"])
name = dag.current_job.name.split("_")
description = dag.current_job.description.split(" ")[:-2] + [f"{natoms} supercell_matrix"]

next_sc = generate_phonopy_job(ctx.settings)
next_sc.name = "_".join(name[:-1] + [str(natoms)])
next_sc.workflow = dag.current_job.workflow
next_sc.description = " ".join(description)

next_sc.wall_time_minutes = update_job["expected_walltime"] / 60

next_sc_data = decode(next_sc.data)
next_sc_data.update(update_data)
next_sc.data = encode(next_sc_data)

next_sc.save()

for child in dag.children:
    dag.add_dependency(next_sc, child)

dag.current_job.children =  []
dag.current_job.save()

dag.add_dependency(dag.current_job, next_sc)
