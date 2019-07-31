"""Functions used to generate a FireWorks Workflow"""

import balsam.launcher.dag as dag

from hilde.balsam.workflow.job_generator import generate_aims_job, generate_phonopy_job

from hilde.helpers.hash import hash_atoms_and_calc
from hilde.phonopy import defaults as ph_defaults


def generate_workflow(workflow_settings):
    """Generates a workflow from given set of steps

    Parameters
    ----------
    workflow_settings: Settings
        The settings object for the desired workflow

    Raises
    ------
    IOError
        If supercell_matrix is not provided for phonopy, phono3py, or statistical_sampling
    """
    jobs = []
    job_dep = {}

    atoms = workflow_settings.atoms.copy()
    atoms.set_calculator(workflow_settings.atoms.calc)

    name = atoms.symbols.get_chemical_formula()
    name += "_" + hash_atoms_and_calc(atoms)[0]

    workflow_settings.general.workdir_local += (
        atoms.symbols.get_chemical_formula() + "/" + hash_atoms_and_calc(atoms)[0] + "/"
    )
    workflow_settings.general.workdir_cluster += (
        atoms.symbols.get_chemical_formula() + "/" + hash_atoms_and_calc(atoms)[0] + "/"
    )

    workflow_settings.atoms = atoms
    if "basisset" in workflow_settings.general:
        basis = workflow_settings.general.pop("basisset")
    else:
        basis = "light"

    # Relaxation
    if workflow_settings.general.get("relax_structure", True):
        jobs.append(
            generate_aims_job(workflow_settings, "light")
        )
        # Tighter Basis Set Relaxation
        use_tight_relax = workflow_settings.general.get("use_tight_relax", False)
        if basis != "light" or use_tight_relax:
            if use_tight_relax:
                basisset_type = "tight"
            else:
                basisset_type = basis
            jobs.append(
                generate_aims_job(workflow_settings, basisset_type)
            )
    for job in jobs:
        job.save()

    for ii in range(len(jobs) - 1):
        job_dep[jobs[ii]] = jobs[ii + 1]

    if len(jobs) > 0:
        final_initialize_job = jobs[-1]
        job_dep[final_initialize_job] = []
    else:
        final_initialize_job = None

    # Phonon Calculations
    phonon_jobs = []
    if "phonopy" in workflow_settings:
        ignore_keys = ["trigonal", "q_mesh"]
        for key, val in ph_defaults.items():
            if key not in workflow_settings.phonopy and key not in ignore_keys:
                workflow_settings.phonopy[key] = val

        workflow_settings.phonopy["basisset_type"] = basis

        if "supercell_matrix" not in workflow_settings.phonopy:
            raise IOError("Initial supercell_matrix must be provided")

        phonon_jobs.append(
            generate_phonopy_job(workflow_settings)
        )
        if final_initialize_job:
            job_dep[final_initialize_job].append(phonon_jobs[0])

    for job in phonon_jobs:
        job.save()

    for key, val in job_dep.items():
        if isinstance(val, list):
            for job in val:
                dag.add_dependency(key, job)
        else:
            dag.add_dependency(key, val)
