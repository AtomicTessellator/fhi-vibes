''' Generate a phonon workflow for Fireworks '''
from fireworks import Workflow
from hilde.fireworks.launchpad import LaunchPadHilde
from hilde.fireworks.workflow_generator import (
    generate_firework,
    get_phonon_serial_task,
    get_phonon_parallel_task,
    get_phonon_analysis_task,
    get_relax_task,
    get_aims_relax_task,
    get_kgrid_task,
)
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.phonopy.wrapper import defaults as ph_defaults
from hilde.phono3py.wrapper import defaults as ph3_defaults

def update_fw_settings(fw_settings, fw_name, queueadapter=None, update_in_spec=True):
    if "out_spec_atoms" in fw_settings and update_in_spec:
        fw_settings["in_spec_atoms"] = fw_settings["out_spec_atoms"]
        fw_settings["in_spec_calc"] = fw_settings["out_spec_calc"]
        fw_settings["from_db"] = True

    fw_settings["out_spec_atoms"] = fw_name + "_atoms"
    fw_settings["out_spec_calc"] = fw_name + "_calc"
    fw_settings["fw_name"] = fw_name
    if "spec" not in fw_settings:
        fw_settings["spec"] = {}
    if queueadapter:
        fw_settings["spec"]["_queueadapter"] = queueadapter
    elif "_queueadapter" in fw_settings["spec"]:
        del(fw_settings["spec"]["_queueadapter"])

    return fw_settings

def generate_fw(atoms, task_list, fw_settings, qadapter, update_settings=None, update_in_spec=True):
    if qadapter:
        fw_settings = update_fw_settings(fw_settings, fw_settings["fw_name"], qadapter, update_in_spec=update_in_spec)
    else:
        fw_settings = update_fw_settings(fw_settings, fw_settings["fw_name"], update_in_spec=update_in_spec)

    fw_settings[
        "fw_name"
    ] += f"_{atoms.symbols.get_chemical_formula()}_{hash_atoms_and_calc(atoms)[0]}"
    if not update_settings:
        update_settings = {}
    at = atoms if "in_spec_atoms" not in  fw_settings else fw_settings["in_spec_atoms"]
    cl = atoms.calc if "in_spec_calc" not in  fw_settings else fw_settings["in_spec_calc"]
    return generate_firework(
        task_list,
        at,
        cl,
        fw_settings.copy(),
        update_calc_settings=update_settings,
    )

def generate_kgrid_fw(atoms, wd, qadapter, fw_settings):
    func_kwargs = {
        "workdir": wd + "/" + fw_settings["fw_name"] + "/",
        "trajectory": "kpt_trajectory.yaml",
    }
    task_spec = get_kgrid_task(func_kwargs)
    return generate_fw(atoms, task_spec, fw_settings, qadapter)

def generate_relax_fw(atoms, wd, fw_settings, qadapter, rel_settings):
    fw_settings["fw_name"] = rel_settings["basisset_type"] + "_relax"
    func_kwargs = {"workdir": wd + "/" + fw_settings["fw_name"] + "/"}
    fw_out_kwargs = {"relax_step": 0}
    task_spec = get_aims_relax_task(func_kwargs, fw_out_kwargs)

    if "rel_method" in rel_settings:
        method = func_kwargs.pop('rel_method')
    else:
        method = "ltrm"

    if "conv_cirt" in rel_settings:
        force_crit = str(func_kwargs.pop("conv_cirt"))
    else:
        force_crit = "5e-3"

    update_settings = {
        "relax_geometry": method + " " + force_crit,
        "relax_unit_cell": "full",
        "basisset_type": rel_settings["basisset_type"],
        "scaled": True,
        "use_sym": True
    }
    return generate_fw(atoms, task_spec, fw_settings, qadapter, update_settings)

def generate_phonon_fw(atoms, wd, fw_settings, qadapter, ph_settings, db_kwargs, update_in_spec=True):
    fw_settings["fw_name"] = ph_settings.pop("type")
    serial = ph_settings.pop("serial")
    update_settings = {}
    if "basisset_type" in ph_settings:
        update_settings["basisset_type"] = ph_settings.pop("basisset_type")
    if "socket_io_port" in ph_settings:
        update_settings["use_pimd_wrapper"] = ph_settings.pop("socket_io_port")
    elif "use_pimd_wrapper" in ph_settings:
        update_settings["use_pimd_wrapper"] = ph_settings.pop("use_pimd_wrapper")
    if "converge_sc" in ph_settings and ph_settings["converge_sc"]:
        func_out_kwargs = dict(db_kwargs, converge_sc=ph_settings.pop("converge_sc"))
    else:
        func_out_kwargs = dict(db_kwargs, converge_sc=False)

    func_kwargs = ph_settings.copy()
    func_kwargs["workdir"] = wd + "/" + fw_settings["fw_name"] + "/"
    if serial:
        task_spec = get_phonon_serial_task(
            "hilde." + fw_settings["fw_name"] + ".workflow.run",
            func_kwargs,
            func_out_kwargs,
            False,
        )
    else:
        task_spec = get_phonon_parallel_task(
            "hilde.tasks.fireworks.preprocess_ph_ph3",
            func_kwargs,
            False,
        )
    return generate_fw(atoms, task_spec, fw_settings, qadapter, update_settings, update_in_spec)

def generate_phonon_postprocess_fw(atoms, wd, fw_settings, ph_settings, db_kwargs):
    if ph_settings['type'] is phonopy:
        fw_settings["mod_spec_add"] = "ph"
    else:
        fw_settings["mod_spec_add"] = "ph3"
    fw_settings["mod_spec_add"] += "_forces"
    fw_settings["fw_name"] = ph_settings.pop("type") + "_analysis"
    serial = ph_settings.pop("serial")

    if "converge_sc" in ph_settings and ph_settings["converge_sc"]:
        func_out_kwargs = dict(db_kwargs, converge_sc=ph_settings.pop("converge_sc"))
    else:
        func_out_kwargs = dict(db_kwargs, converge_sc=False)

    func_kwargs = ph_settings.copy()
    func_kwargs["workdir"] = wd + "/" + fw_settings["fw_name"] + "/"

    task_spec = get_phonon_analysis_task(
        "hilde."+ fw_settings["fw_name"][:-9] + ".postprocess.postprocess",
        func_kwargs,
        fw_settings,
        "metadata_"+fw_settings["mod_spec_add"][:-7],
        func_out_kwargs,
    )
    return generate_fw(atoms, task_spec, fw_settings, None, update_settings)

def generate_phonon_workflow(workflow, atoms, fw_settings):
    '''
    Generates a workflow from given set of steps
    Args
        steps (list of dicts): List of parameters for all the steps in a given system
        atoms (ASE atoms object, dict): ASE Atoms object to preform the calculation on
    Returns (Workflow or None):
        Either adds the workflow to the launchpad or returns it
    '''
    fw_steps = []
    fw_dep = {}

    fw_settings["fw_name"] = "kgrid_opt"
    fw_settings["out_spec_k_den"] = "kgrid"
    # K-grid optimization
    if "kgrid_qadapter" in workflow:
        qadapter = workflow["kgrid_qadapter"]
    else:
        qadapter = None
    fw_steps.append(
        generate_kgrid_fw(
            atoms,
            workflow.general.workdir_cluster,
            qadapter,
            fw_settings,
        )
    )
    # Light Basis Set Relaxation
    fw_settings["kpoint_density_spec"] = "kgrid"
    del(fw_settings["out_spec_k_den"])
    light_relax_set = { "basisset_type": "light"}
    if "light_rel_qadapter" in workflow:
        qadapter = workflow["light_rel_qadapter"]
    else:
        qadapter = None

    fw_steps.append(
        generate_relax_fw(
            atoms,
            workflow.general.workdir_cluster,
            fw_settings,
            qadapter,
            light_relax_set,
        )
    )
    fw_dep[fw_steps[-2]] = fw_steps[-1]

    # Tight Basis Set Relaxation
    tight_relax_set = { "basisset_type": "tight"}
    if "tight_rel_qadapter" in workflow:
        qadapter = workflow["tight_rel_qadapter"]
    else:
        qadapter = None
    fw_steps.append(
        generate_relax_fw(
            atoms,
            workflow.general.workdir_cluster,
            fw_settings,
            qadapter,
            tight_relax_set
        )
    )
    fw_dep[fw_steps[-2]] = fw_steps[-1]

    pre_ph_fw = fw_steps[-1]
    fw_dep[pre_ph_fw] = []

    # Phonon Calculations
    if "db_storage" in workflow:
        db_kwargs = workflow.db_storage
    else:
        db_kwargs = None
    phonopy_set = ph_defaults.copy()
    del(phonopy_set['trigonal'])
    del(phonopy_set['q_mesh'])
    phonopy_set ['serial'] = True
    phonopy_set ['type'] = "phonopy"

    if "phonopy_qadapter" in workflow:
        qadapter = workflow["phonopy_qadapter"]
    else:
        qadapter = None
    if "phonopy" in workflow:
        for key,val in workflow.phonopy.items():
            phonopy_set[key] = val
    if "supercell_matrix" not in phonopy_set:
        phonopy_set["supercell_matrix"] = np.eye(3)
        phonopy_set["converge_sc"] = True
    fw_steps.append(
        generate_phonon_fw(
            atoms,
            workflow.general.workdir_cluster,
            fw_settings,
            qadapter,
            phonopy_set.copy(),
            db_kwargs,
        )
    )
    fw_dep[pre_ph_fw].append(fw_steps[-1])
    if not phonopy_set["serial"]:
        fw_steps.append(
            generate_phonon_postprocess_fw(
                atoms,
                workflow.workdir_local,
                fw_settings,
                phonopy_set,
                db_kwargs,
            )
        )
        fw_dep[fw_steps[-2]] = fw_steps[-1]

    if ("third_order" in workflow.general and workflow.general.third_order) or "phono3py" in workflow:
        phono3py_set = ph3_defaults.copy()
        phono3py_set['serial'] = True
        phono3py_set ['type'] = "phono3py"
        del(phono3py_set["displacement"])
        del(phono3py_set["cutoff_pair_distance"])
        del(phono3py_set["q_mesh"])
        if "phono3py_qadapter" in workflow:
            qadapter = workflow["phono3py_qadapter"]
        else:
            qadapter = None
        if "phono3py" in workflow:
            for key,val in workflow.phono3py.items():
                phono3py_set[key] = val
        if "supercell_matrix" not in phono3py_set:
            phono3py_set["supercell_matrix"] = np.eye(3)
            phono3py_set["converge_sc"] = True
        fw_steps.append(
            generate_phonon_fw(
                atoms,
                workflow.general.workdir_cluster,
                fw_settings,
                qadapter,
                phono3py_set.copy(),
                db_kwargs,
                update_in_spec=False,
            )
        )
        fw_dep[pre_ph_fw].append(fw_steps[-1])
        if not phono3py_set["serial"]:
            fw_steps.append(
                generate_phonon_postprocess_fw(
                    atoms,
                    workflow.workdir_local,
                    fw_settings,
                    workflow.phono3py.copy(),
                    db_kwargs,
                )
            )
            fw_dep[fw_steps[-2]] = fw_steps[-1]
    if "launchpad_yaml" in fw_settings:
        launchpad = LaunchPadHilde.from_file(fw_settings["launchpad_yaml"])
    else:
        launchpad = LaunchPadHilde.auto_load()
    launchpad.add_wf(Workflow(fw_steps, fw_dep, name=fw_settings["name"]))
    return None
