"""Functions used to setup phonopy postprocessing"""

from glob import glob

from pathlib import Path

from ase.io.aims import read_aims_output

import balsam.launcher.dag as dag

from hilde.aims.context import AimsContext
from hilde.aims.setup import setup_aims

from hilde.balsam.workflow.job_generator import (
    generate_gruneisen_jobs,
    generate_stat_samp_jobs,
)
from hilde.phonopy import displacement_id_str
from hilde.phonopy.context import PhonopyContext
from hilde.settings import Settings
from hilde.trajectory import metadata2file, step2file


def collect_calcs_to_trajectory(calc_direc, metadata):
    """Collect a series of calculations into a trajectory.son file"""
    trajectory = Path(calc_direc) / "trajectory.son"

    files = list(glob(f"{calc_direc}/*/aims.out"))
    files = sorted(files, key=lambda s: int(s.split("/")[-2]))
    calculated_atoms = [read_aims_output(file) for file in files]
    # Add the displacement ids
    for file, atoms in zip(files, calculated_atoms):
        if atoms.info is None:
            atoms.info = dict()
        atoms.info[displacement_id_str] = int(file.split("/")[-2])

    # Create the trajectory file
    metadata2file(metadata, trajectory)
    for atoms in calculated_atoms:
        step2file(atoms, atoms.calc, trajectory)


def setup_gruneisen(vol_factor, sc_mat, walltime):
    """Setup a second structure for a gruneisen FD calculation"""
    ctx = PhonopyContext(
        Settings(settings_file="phonopy.in", read_config=False), read_config=False
    )
    calc = setup_aims(AimsContext(Settings(settings_file="phonopy.in")))
    calc.parameters.pop("use_pimd_wrapper", None)

    ctx.settings.atoms.set_calculator(calc)

    ctx.settings.phonopy["supercell_matrix"] = sc_mat
    ctx.settings.obj["supercell_matrix"] = sc_mat

    ctx.settings.phonopy["get_gruniesen"] = False
    ctx.settings.phonopy["converge_phonons"] = False
    ctx.settings.pop("statistical_sampling", None)

    ctx.settings.write("gruneisen.in")
    ctx.settings._settings_file = "gruneisen.in"

    jobs = generate_gruneisen_jobs(ctx.settings, vol_factor)

    Path("gruneisen.in").unlink()

    for job in jobs:
        dag.add_dependency(dag.current_job, job)
        job.wall_time_minues = walltime
        job.save()


def setup_statistical_sampling(sc_mat, phonon, walltime):
    """Setup a statistical sampling run"""
    ctx = PhonopyContext(
        Settings(settings_file="phonopy.in", read_config=False), read_config=False
    )
    settings = ctx.settings
    calc = setup_aims(AimsContext(Settings(settings_file="phonopy.in")))
    calc.parameters.pop("use_pimd_wrapper", None)

    settings.atoms.set_calculator(calc)
    settings.statistical_sampling["supercell_matrix"] = sc_mat

    settings.general["relax_structure"] = False
    settings.pop("relaxation", None)

    ctx.settings.write("stat_samp.in")
    ctx.settings._settings_file = "stat_samp.in"

    jobs = generate_stat_samp_jobs(settings, phonon)

    Path("stat_samp.in").unlink()

    for job in jobs:
        dag.add_dependency(dag.current_job, job)
        job.wall_time_minues = walltime
        job.save()
