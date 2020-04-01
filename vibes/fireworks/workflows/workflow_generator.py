"""Functions used to generate a FireWorks Workflow"""
from pathlib import Path

import numpy as np
from fireworks import Workflow
from vibes.fireworks.launchpad import LaunchPad
from vibes.fireworks.workflows.firework_generator import (
    generate_aims_fw,
    generate_aims_relax_fw,
    generate_gruniesen_fd_fw,
    generate_kgrid_fw,
    generate_md_fw,
    generate_phonon_fw,
    generate_phonon_postprocess_fw,
    generate_relax_fw,
    generate_stat_samp_fw,
    generate_stat_samp_postprocess_fw,
)
from vibes.helpers.hash import hash_atoms_and_calc
from vibes.helpers.numerics import get_3x3_matrix
from vibes.phonopy._defaults import kwargs as ph_defaults


def process_aims_relaxation(workflow, atoms, fw_settings, basis):
    """Processes the workflow settings to get all relaxation steps

    Parameters
    ----------
    workflow : settings.Settings
        Settings for the workflow
    atoms : ase.atoms.Atoms
        The Atoms object for the structure
    fw_settings : dict
        FireWorks specific settings
    basis : str
        The default basis set used for this calculation:

    Returns
    -------
    list of FireWorks
        Relaxation FireWorks to add to the workflow

    """
    fw_steps = []
    fw_steps.append(
        generate_aims_relax_fw(workflow.settings, atoms, fw_settings, "light")
    )

    # Tighter Basis Set Relaxation
    use_tight_relax = workflow.settings.general.get("use_tight_relax", False)

    if basis != "light" or use_tight_relax:
        if use_tight_relax:
            basisset_type = "tight"
        else:
            basisset_type = basis
        fw_steps.append(
            generate_aims_relax_fw(workflow.settings, atoms, fw_settings, basisset_type)
        )

    return fw_steps


def process_phonons(workflow, atoms, fw_settings, basis):
    """Processes the workflow settings to get all phonopy steps

    Parameters
    ----------
    workflow : settings.Settings
        Settings for the workflow
    atoms : ase.atoms.Atoms
        The Atoms object for the structure
    fw_settings : dict
        FireWorks specific settings
    basis : str
        The default basis set used for this calculation

    Returns
    -------
    list of FireWorks
        Phonopy related FireWorks to add to the workflow

    Raises
    ------
    ValueError
        If supercell_matrix is not provided for phonopy, phono3py, or statistical_sampling

    """
    phonon_fws = []
    ignore_keys = ["trigonal", "q_mesh"]
    for key, val in ph_defaults.items():
        if key not in workflow.settings.phonopy and key not in ignore_keys:
            workflow.settings.phonopy[key] = val

    if "serial" not in workflow.settings.phonopy:
        workflow.settings.phonopy["serial"] = True

    if basis:
        workflow.settings.phonopy["basisset_type"] = basis
    else:
        workflow.settings.phonopy.pop("basisset_type", None)

    if "supercell_matrix" not in workflow.settings.phonopy:
        raise ValueError("Initial supercell_matrix must be provided")

    phonon_fws.append(
        generate_phonon_fw(workflow.settings, atoms, fw_settings, "phonopy")
    )
    phonon_fws.append(
        generate_phonon_postprocess_fw(workflow.settings, atoms, fw_settings, "phonopy")
    )
    if getattr(workflow.settings.phonopy, "get_gruniesen", False):
        phonon_fws += process_grun(workflow, atoms, fw_settings)
    return phonon_fws


def process_stat_samp(workflow, atoms, fw_settings):
    """Processes the workflow settings to get all Statistical Sampling steps

    Parameters
    ----------
    workflow : settings.Settings
        Settings for the workflow
    atoms : ase.atoms.Atoms
        The Atoms object for the structure
    fw_settings : dict
        FireWorks specific settings

    Returns
    -------
    list of FireWorks
        Statistical Sampling related FireWorks to add to the workflow

    """
    stat_samp_fws = []
    stat_samp_fws.append(generate_stat_samp_fw(workflow.settings, atoms, fw_settings))

    stat_samp_fws.append(
        generate_stat_samp_postprocess_fw(workflow.settings, atoms, fw_settings)
    )
    return stat_samp_fws


def process_grun(workflow, atoms, fw_settings):
    """Processes the workflow settings to get all Gruneisen steps

    Parameters
    ----------
    workflow : settings.Settings
        Settings for the workflow
    atoms : ase.atoms.Atoms
        The Atoms object for the structure
    fw_settings : dict
        FireWorks specific settings

    Returns
    -------
    list of FireWorks
        Gruneisen parameter related FireWorks to add to the workflow

    """
    grun_fws = []
    if getattr(workflow.settings.phonopy, "converge_phonons", False):
        trajectory_file = (
            workflow.settings.general.workdir_local + "/converged/trajectory.son"
        )
    else:
        sc_mat = get_3x3_matrix(workflow.settings.phonopy.supercell_matrix)
        natoms = int(round(np.linalg.det(sc_mat) * len(atoms)))
        trajectory_file = (
            workflow.settings.general.workdir_local
            + f"/sc_natoms_{natoms}/phonopy_analysis/trajectory.son"
        )

    constraints = []
    for cosntr in atoms.constraints:
        constraints.append(cosntr.todict())

    grun_fws.append(
        generate_gruniesen_fd_fw(
            workflow.settings, atoms, trajectory_file, constraints, fw_settings
        )
    )
    return grun_fws


def process_workdir(workdir, atoms, make_absolute):
    """process the working directory

    Parameters
    ----------
    workdir: str
        input working directory
    atoms: ase.atoms.Atoms
        ASE Atoms object for the working directory
    make_absolute: bool
        If True make the working directory absolute

    Returns
    -------
    workdir: str
        output working directory

    """
    if make_absolute:
        workdir = Path(workdir).absolute()

    workdir = (
        workdir
        / atoms.symbols.get_chemical_formula(mode="metal", empirical=True)
        / hash_atoms_and_calc(atoms)[0]
    )

    return str(workdir) + "/"


def generate_workflow(workflow, atoms, launchpad_yaml=None, make_absolute=True):
    """Generates a workflow from given set of steps

    Parameters
    ----------
    workflow : Settings
        The settings object for the desired workflow
    atoms : ase.atoms.Atoms
        ASE Atoms object to preform the calculation on, with an attached calculator
    launchpad_yaml : str
        filename for the launchpad definition file (Default value = None)
    make_absolute: bool
        If True make the working directory absolute

    Returns
    -------
    fireworks.Workflow
        The FireWorks Workflow for the given workflow

    """
    fw_steps = []
    fw_dep = {}

    fw_settings = {
        "name": atoms.symbols.get_chemical_formula(mode="metal", empirical=True)
        + "_"
        + hash_atoms_and_calc(atoms)[0]
    }
    wd_cluster = workflow.settings.general.get(
        "workdir_cluster", workflow.settings.general.workdir_local
    )

    workflow.settings.general["workdir_local"] = process_workdir(
        workflow.settings.general.workdir_local, atoms, make_absolute
    )
    workflow.settings.general["workdir_cluster"] = process_workdir(
        wd_cluster, atoms, make_absolute
    )
    if atoms.calc.name == "aims":
        # K-grid optimization
        if workflow.settings.general.get("opt_kgrid", False):
            fw_steps.append(generate_kgrid_fw(workflow.settings, atoms, fw_settings))

        basis = workflow.settings.general.get("basisset", "light")

        # Relaxation
        if workflow.settings.general.get("relax_structure", False):
            fw_steps += process_aims_relaxation(workflow, atoms, fw_settings, basis)
    else:
        basis = None
        if "relaxation" in workflow.settings:
            fw_steps.append(generate_relax_fw(workflow.settings, atoms, fw_settings))

    # Setup workflow branching point
    for ii in range(len(fw_steps) - 1):
        fw_dep[fw_steps[ii]] = fw_steps[ii + 1]

    if fw_steps:
        final_initialize_fw = fw_steps[-1]
        fw_dep[final_initialize_fw] = []
    else:
        final_initialize_fw = None

    # Phonon Calculations
    if "phonopy" in workflow.settings:
        phonon_fws = process_phonons(workflow, atoms, fw_settings, basis)
        fw_dep[phonon_fws[0]] = phonon_fws[1]
        if final_initialize_fw:
            fw_dep[final_initialize_fw].append(phonon_fws[0])

        if len(phonon_fws) > 2:
            fw_dep[phonon_fws[1]] = [phonon_fws[2]]
        else:
            fw_dep[phonon_fws[1]] = []
        fw_steps += phonon_fws

    # Statistical Sampling
    if "statistical_sampling" in workflow.settings:
        stat_samp_fws = process_stat_samp(workflow, atoms, fw_settings)

        if "phonopy" in workflow.settings:
            fw_dep[phonon_fws[1]].append(stat_samp_fws[0])
        elif final_initialize_fw:
            fw_dep[final_initialize_fw].append(stat_samp_fws[0])

        fw_dep[stat_samp_fws[0]] = stat_samp_fws[1]
        fw_steps += stat_samp_fws

    # Molecular dynamics
    if "md" in workflow.settings:
        md_fws = generate_md_fw(workflow.settings, atoms, fw_settings)
        if "phonopy" in workflow.settings:
            fw_dep[phonon_fws[1]] += md_fws
        elif final_initialize_fw:
            fw_dep[final_initialize_fw] += md_fws
        fw_steps += md_fws

    # Aims Calculations if no other term is present
    if not fw_steps:
        fw_steps.append(generate_aims_fw(workflow.settings, atoms, fw_settings))

    if launchpad_yaml:
        launchpad = LaunchPad.from_file(launchpad_yaml)
        launchpad.add_wf(Workflow(fw_steps, fw_dep, name=fw_settings["name"]))
        return None

    return Workflow(fw_steps, fw_dep, name=fw_settings["name"])
