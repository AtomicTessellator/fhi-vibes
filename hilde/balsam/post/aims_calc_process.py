#!/usr/bin/env python
"""Check the completion of an aims calculation"""
from pathlib import Path
from glob import glob

import balsam.launcher.dag as dag

from hilde.aims.context import AimsContext
from hilde.aims.setup import setup_aims

from hilde.balsam.data_encoder import decode, encode
from hilde.balsam.workflow.job_generator import generate_aims_job
from hilde.helpers.attribute_dict import AttributeDict
from hilde.helpers.k_grid import k2d
from hilde.phonon_db.ase_converters import calc2dict, dict2atoms
from hilde.settings import Settings
from hilde.fireworks.tasks.postprocess.relax import check_aims


def postprocess_aims(data):
    Path(glob("aims_calc_*")[0]).rename("aims.out")

    workdir = "./"

    ctx = AimsContext(Settings(settings_file="aims.in"))
    calc = setup_aims(ctx)

    atoms = ctx.settings.get_atoms()
    atoms.set_calculator(calc)
    atoms.info.update(data["atoms"]["info"])

    results = check_aims(
        data["atoms"],
        calc,
        atoms,
        calc_number=data["calc_number"],
        walltime=dag.current_job.wall_time_minutes * 60,
        workdir=workdir,
    )
    completed, calc_number, new_atoms_dict, walltime = results

    if not completed and walltime > 0:
        calc.parameters["walltime"] = walltime

    calc_dict = calc2dict(calc)

    k_pt_density = k2d(
        dict2atoms(new_atoms_dict), calc_dict["calculator_parameters"]["k_grid"]
    )

    calc_dict.pop("results", None)
    calc_dict["calculator_parameters"].pop("k_grid", None)
    calc_dict["calculator_parameters"].pop("use_pimd_wrapper", None)

    if completed:
        calc_dict["calculator_parameters"].pop("relax_geometry", None)
        calc_dict["calculator_parameters"].pop("relax_unit_cell", None)

    update_data = {
        "atoms": new_atoms_dict,
        "calculator": calc_dict,
        "k_pt_density": k_pt_density,
    }

    data.update(update_data)
    for child in dag.children:
        ch_data = decode(child.data)
        ch_data.update(update_data)
        child.data = encode(ch_data.copy())
        child.save()

    if completed:
        return True, data

    settings = AimsContext(
        Settings(settings_file="aims.in", read_config=False)
    ).settings
    settings.atoms = dict2atoms(new_atoms_dict)
    settings["control"].update(calc_dict["calculator_parameters"])
    settings["control_kpt"] = AttributeDict()
    settings["control_kpt"]["density"] = k_pt_density

    name = dag.current_job.name.split("_")

    restart = generate_aims_job(settings, calc_number=calc_number)
    restart.name = "_".join(name[:-1] + [calc_number])
    restart.workflow = dag.current_job.workflow
    restart.description = (
        dag.current_job.description[:-1] + f" restart number {calc_number}"
    )
    restart.application = dag.current_job.application
    restart_data = decode(restart.data)
    for key, val in data.items():
        if key not in restart_data:
            restart_data[key] = val
    restart.data = encode(restart_data)
    if walltime > 0:
        restart.wall_time_minutes = walltime / 60

    dag.add_dependency(dag.current_job, restart)

    for child in dag.children:
        dag.add_dependency(restart, child)

    restart.save()

    return False, data
