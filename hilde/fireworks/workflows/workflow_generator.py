"""Functions used to generate a FireWorks Workflow"""
from fireworks import Workflow

from hilde.fireworks.launchpad import LaunchPad
from hilde.fireworks.workflows.firework_generator import (
    generate_kgrid_fw,
    generate_relax_fw,
    generate_phonon_fw,
    generate_phonon_postprocess_fw,
    # generate_stat_samp_fw,
    generate_aims_fw,
)
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.phonopy import defaults as ph_defaults


def generate_workflow(workflow_settings, atoms, launchpad_yaml=None):
    """Generates a workflow from given set of steps

    Parameters
    ----------
    workflow_settings: Settings
        The settings object for the desired workflow
    atoms: ase.atoms.Atoms
        ASE Atoms object to preform the calculation on, with an attached calculator
    launchpad_yaml: str
        filename for the launchpad definition file

    Raises
    ------
    IOError
        If supercell_matrix is not provided for phonopy, phono3py, or statistical_sampling
    """
    fw_steps = []
    fw_dep = {}

    fw_settings = {
        "name": atoms.symbols.get_chemical_formula()
        + "_"
        + hash_atoms_and_calc(atoms)[0]
    }
    workflow_settings.general.workdir_local += (
        atoms.symbols.get_chemical_formula() + "/" + hash_atoms_and_calc(atoms)[0] + "/"
    )
    workflow_settings.general.workdir_cluster += (
        atoms.symbols.get_chemical_formula() + "/" + hash_atoms_and_calc(atoms)[0] + "/"
    )

    # K-grid optimization
    if workflow_settings.general.get("opt_kgrid", True):
        fw_steps.append(generate_kgrid_fw(workflow_settings, atoms, fw_settings))

        # Now that the k-grid has been optimized use that data
        fw_settings["kpoint_density_spec"] = "kgrid"
        del fw_settings["out_spec_k_den"]

    if "basisset" in workflow_settings.general:
        basis = workflow_settings.general.pop("basisset")
    else:
        basis = "light"

    # Relaxation
    if workflow_settings.general.get("relax_structure", True):
        fw_steps.append(
            generate_relax_fw(workflow_settings, atoms, fw_settings, "light")
        )
        # Tighter Basis Set Relaxation
        use_tight_relax = workflow_settings.general.get("use_tight_relax", False)
        if basis != "light" or use_tight_relax:
            if use_tight_relax:
                basisset_type = "tight"
            else:
                basisset_type = basis
            fw_steps.append(
                generate_relax_fw(workflow_settings, atoms, fw_settings, basisset_type)
            )
    for ii in range(len(fw_steps) - 1):
        fw_dep[fw_steps[ii]] = fw_steps[ii + 1]

    if len(fw_steps) > 0:
        final_initialize_fw = fw_steps[-1]
        fw_dep[final_initialize_fw] = []
    else:
        final_initialize_fw = None

    # Phonon Calculations
    phonon_fws = []
    phonon3_fws = []
    stat_samp_fws = []
    aims_calc_fws = []
    if "phonopy" in workflow_settings:
        ignore_keys = ["trigonal", "q_mesh"]
        for key, val in ph_defaults.items():
            if key not in workflow_settings.phonopy and key not in ignore_keys:
                workflow_settings.phonopy[key] = val

        if "serial" not in workflow_settings.phonopy:
            workflow_settings.phonopy["serial"] = True

        workflow_settings.phonopy["basisset_type"] = basis

        if "supercell_matrix" not in workflow_settings.phonopy:
            raise IOError("Initial supercell_matrix must be provided")

        phonon_fws.append(
            generate_phonon_fw(workflow_settings, atoms, fw_settings, "phonopy")
        )
        phonon_fws.append(
            generate_phonon_postprocess_fw(
                workflow_settings, atoms, fw_settings, "phonopy"
            )
        )
        if final_initialize_fw:
            fw_dep[final_initialize_fw].append(phonon_fws[0])
        fw_dep[phonon_fws[0]] = phonon_fws[1]

    if "phono3py" in workflow_settings:
        from hilde.phono3py import defaults as ph3_defaults

        ignore_keys = ["displacement", "cutoff_pair_distance", "q_mesh"]
        for key, val in ph3_defaults.items():
            if key not in workflow_settings.phono3py and key not in ignore_keys:
                workflow_settings.phono3py[key] = val

        if "serial" not in workflow_settings.phonopy:
            workflow_settings.phono3py["serial"] = True

        workflow_settings.phono3py["basisset_type"] = basis

        if "supercell_matrix" not in workflow_settings.phono3py:
            raise IOError("Initial supercell_matrix must be provided")

        phonon3_fws.append(
            generate_phonon_fw(workflow_settings, atoms, fw_settings, "phono3py")
        )

        phonon3_fws.append(
            generate_phonon_postprocess_fw(
                workflow_settings, atoms, fw_settings, "phono3py"
            )
        )
        if final_initialize_fw:
            fw_dep[final_initialize_fw].append(phonon3_fws[0])
        fw_dep[phonon3_fws[0]] = phonon3_fws[1]

    # Harmonic Analysis
    # if "statistical_sampling" in workflow_settings:
    #     if "phonopy" in workflow_settings:
    #         if workflow_settings.phonopy.get("converge_phonons", False):
    #             workflow_settings.statistical_sampling[
    #                 "phonon_file"
    #             ] = f"{workflow_settings.general.workdir_local}/converged/trajectory.son"
    #         else:
    #             workflow_settings.statistical_sampling[
    #                 "phonon_file"
    #             ] = f"{workflow_settings.general.workdir_local}//phonopy_analysis/trajectory.son"
    #     stat_samp_fws.append(
    #         generate_stat_samp_fw(workflow_settings, atoms, fw_settings)
    #     )
    #     if "phonopy" in workflow_settings:
    #         fw_dep[phonon_fws[-1]] = stat_samp_fws[0]
    #     elif final_initialize_fw:
    #         fw_dep[final_initialize_fw].append(stat_samp_fws[0])

    # Aims Calculations if no other term is present
    if not fw_steps and not phonon_fws and not phonon3_fws and not stat_samp_fws:
        aims_calc_fws.append(generate_aims_fw(workflow_settings, atoms, fw_settings))

    for fw in phonon_fws:
        fw_steps.append(fw)
    for fw in phonon3_fws:
        fw_steps.append(fw)
    for fw in stat_samp_fws:
        fw_steps.append(fw)
    for fw in aims_calc_fws:
        fw_steps.append(fw)

    if launchpad_yaml:
        launchpad = LaunchPad.from_file(launchpad_yaml)
    else:
        launchpad = LaunchPad.auto_load()
    launchpad.add_wf(Workflow(fw_steps, fw_dep, name=fw_settings["name"]))
