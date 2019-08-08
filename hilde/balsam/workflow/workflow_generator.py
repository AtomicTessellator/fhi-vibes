"""Functions used to generate a FireWorks Workflow"""

from pathlib import Path

import balsam.launcher.dag as dag

from hilde.balsam.data_encoder import encode, decode
from hilde.balsam.workflow.job_generator import (
    generate_aims_job,
    make_job_phonopy_setup,
    generate_phonopy_jobs,
    get_phonon_setup_data,
)
from hilde.helpers.attribute_dict import AttributeDict
from hilde.helpers.hash import hash_atoms
from hilde.helpers.numerics import get_3x3_matrix
from hilde.phonon_db.ase_converters import atoms2dict
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

    atoms = workflow_settings.atoms.copy()
    atoms.set_calculator(workflow_settings.atoms.calc)

    name = atoms.symbols.get_chemical_formula()
    name += "_" + hash_atoms(atoms)[0]

    workflow_settings.atoms = atoms
    if "basisset" in workflow_settings.general:
        basis = workflow_settings.general.pop("basisset")
    else:
        basis = "light"

    # Phonon Calculations setup
    if "phonopy" in workflow_settings:
        if "supercell_matrix" not in workflow_settings.phonopy:
            raise IOError("Initial supercell_matrix must be provided")

        workflow_settings.phonopy["supercell_matrix"] = get_3x3_matrix(
            workflow_settings.phonopy.supercell_matrix
        )

        ignore_keys = ["trigonal", "q_mesh"]
        for key, val in ph_defaults.items():
            if key not in workflow_settings.phonopy and key not in ignore_keys:
                workflow_settings.phonopy[key] = val

        workflow_settings.phonopy["basisset_type"] = basis

    # Relaxation
    if workflow_settings.general.get("relax_structure", True):
        jobs.append(generate_aims_job(workflow_settings, "light"))
        # Tighter Basis Set Relaxation
        use_tight_relax = workflow_settings.general.get("use_tight_relax", False)
        if basis != "light" or use_tight_relax:
            if use_tight_relax:
                basisset_type = "tight"
            else:
                basisset_type = basis
            jobs.append(generate_aims_job(workflow_settings, basisset_type))

    for job in jobs:
        job.save()

    for ii in range(len(jobs) - 1):
        dag.add_dependency(jobs[ii], jobs[ii + 1])

    if "phonopy" in workflow_settings:
        if jobs:
            make_job_phonopy_setup(workflow_settings, jobs[-1])
        else:
            ph_data = get_phonon_setup_data(
                workflow_settings, workflow_settings.pop("phonopy_qadapter")
            )
            sc_matrix = workflow_settings.phonopy.get("supercell_matrix")
            ph_data["sc_matrix_original"] = sc_matrix
            del ph_data["ph_primitive"]
            ph_data["conv_crit"] = workflow_settings.phonopy.get("conv_crit", 0.95)
            ph_data["ctx"] = workflow_settings
            data = dict(ph_data=ph_data)
            unlink_geo = True
            if workflow_settings.geometry.get("file", None) != "geometry.in":
                workflow_settings["geometry"] = AttributeDict({"file": "geometry.in"})
                unlink_geo = False
            atoms.write("geometry.in", scaled=True, geo_constrain=True)

            workflow_settings.pop("restart", None)
            workflow_settings.pop("relaxation", None)
            workflow_settings.general["relax_structure"] = False
            jobs = generate_phonopy_jobs(workflow_settings, data)
            if unlink_geo:
                Path("geometry.in").unlink()
            for job in jobs:
                job_data = decode(job.data)
                job_data["primitive"] = atoms2dict(atoms)
                job.data = encode(job_data)
                job.save()
