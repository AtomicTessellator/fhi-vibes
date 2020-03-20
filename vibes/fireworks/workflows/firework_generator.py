"""Generates Task Specific FireWorks"""
import numpy as np
from fireworks import Firework

from vibes.filenames import filenames
from vibes.fireworks.tasks.postprocess.phonons import time2str
from vibes.fireworks.tasks.task_spec import TaskSpec
from vibes.fireworks.tasks.utility_tasks import update_calc
from vibes.fireworks.workflows.task_generator import (
    generate_mod_calc_task,
    generate_task,
    generate_update_calc_task,
)
from vibes.fireworks.workflows.task_spec_generator import (
    gen_aims_task_spec,
    gen_gruniesen_task_spec,
    gen_kgrid_task_spec,
    gen_md_task_spec,
    gen_phonon_analysis_task_spec,
    gen_phonon_task_spec,
    gen_stat_samp_analysis_task_spec,
    gen_stat_samp_task_spec,
)
from vibes.helpers.converters import atoms2dict, calc2dict
from vibes.helpers.hash import hash_atoms_and_calc
from vibes.helpers.k_grid import k2d, update_k_grid, update_k_grid_calc_dict
from vibes.helpers.numerics import get_3x3_matrix
from vibes.helpers.watchdogs import str2time
from vibes.phonopy.wrapper import preprocess


def update_fw_settings(fw_settings, fw_name, queueadapter=None, update_in_spec=True):
    """update the fw_settings for the next step

    Parameters
    ----------
    fw_settings:(ict
        Current fw_settings
    fw_name:(tr
        name of the current step
    queueadapter:(ict
        dict describing the queueadapter changes for this firework
    update_in_spec:(ool
        If true move current out_spec to be in_spec

    Returns
    -------
    dict
        The updated fw_settings
    """
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
        del fw_settings["spec"]["_queueadapter"]

    return fw_settings


def generate_firework(
    task_spec_list=None,
    atoms=None,
    calculator=None,
    fw_settings=None,
    atoms_calc_from_spec=False,
    update_calc_settings=None,
    func=None,
    func_fw_out=None,
    func_kwargs=None,
    func_fw_out_kwargs=None,
    args=None,
    inputs=None,
):
    """A function that takes in a set of inputs and returns a Firework for them

    Parameters
    ----------
    task_spec_list: list of TaskSpecs
        list of task specifications to perform
    atoms: ase.atoms.Atoms, dictionary or str
        If not atoms_calc_from_spec then this must be an ASE Atoms object or a dict
        If atoms_calc_from_spec then this must be a key str
    calculator: ase.calculators.calulator.Calculator, dictionary or str
        If not atoms_calc_from_spec then this must be an ASE Calculator or a dict
        If atoms_calc_from_spec then this must be a key str
    fw_settings: dict
        Settings used by fireworks to place objects in the right part of the MongoDB
    atoms_calc_from_spec: bool
        If True retrieve the atoms/Calculator objects from the MongoDB launchpad
    update_calc_settings: dict
        Used to update the Calculator parameters
    func: str
        Function path for the firework
    func_fw_out: str
        Function path for the fireworks FWAction generator
    func_kwargs: dict
        Keyword arguments for the main function
    func_fw_out_kwargs: dict
        Keyword arguments for the fw_out function
    args: list
        List of arguments to pass to func
    inputs:(ist
        List of spec to pull in as args from the FireWorks Database

    Returns
    -------
    Firework
        A Firework that will perform the desired operation on a set of atoms

    Raises
    ------
    IOError
        If conflicting task_spec definitions provided, or none are provided
    """
    fw_settings = fw_settings.copy()
    if "spec" not in fw_settings:
        fw_settings["spec"] = {}

    if update_calc_settings is None:
        update_calc_settings = {}

    if func:
        if task_spec_list:
            raise IOError(
                "You have defined both a task_spec_list and arguments to generate one"
            )
        task_with_atoms_obj = atoms is not None
        task_spec_list = [
            TaskSpec(
                func,
                func_fw_out,
                task_with_atoms_obj,
                func_kwargs,
                func_fw_out_kwargs,
                args,
                inputs,
            )
        ]
    elif not task_spec_list:
        raise IOError(
            "You have not defined a task_spec_list or arguments to generate one"
        )
    if isinstance(task_spec_list, TaskSpec):
        task_spec_list = [task_spec_list]

    atoms_calc_from_spec = fw_settings.get("from_db", False)

    if "fw_name" not in fw_settings:
        fw_settings["fw_base_name"] = ""
    elif "fw_base_name" not in fw_settings:
        fw_settings["fw_base_name"] = fw_settings["fw_name"]

    setup_tasks = []
    if atoms:
        if not atoms_calc_from_spec:
            # Preform calculator updates here
            at = atoms2dict(atoms, add_constraints=True)

            if not isinstance(calculator, str):
                if "k_grid_density" in update_calc_settings:
                    if not isinstance(calculator, dict):
                        update_k_grid(
                            atoms, calculator, update_calc_settings["k_grid_density"]
                        )
                    else:
                        recipcell = np.linalg.pinv(at["cell"]).transpose()
                        calculator = update_k_grid_calc_dict(
                            calculator,
                            recipcell,
                            update_calc_settings["k_grid_density"],
                        )

                cl = calc2dict(calculator)

                for key, val in update_calc_settings.items():
                    if key not in ("k_grid_density", "kgrid"):
                        cl = update_calc(cl, key, val)
                if cl["calculator"].lower() == "aims":
                    fw_settings["spec"]["kgrid"] = k2d(
                        atoms, cl["calculator_parameters"]["k_grid"]
                    )
                else:
                    fw_settings["spec"]["kgrid"] = None
            else:
                cl = calculator
                setup_tasks.append(
                    generate_update_calc_task(calculator, update_calc_settings)
                )
        else:
            # Add tasks to update calculator parameters
            at = atoms
            cl = calculator
            if update_calc_settings.keys():
                setup_tasks.append(
                    generate_update_calc_task(calculator, update_calc_settings)
                )

            setup_tasks.append(generate_mod_calc_task(at, cl, "calculator", "kgrid"))
            cl = "calculator"
    else:
        at = None
        cl = None
    job_tasks = [generate_task(ts, fw_settings, at, cl) for ts in task_spec_list]

    return Firework(
        setup_tasks + job_tasks, name=fw_settings["fw_name"], spec=fw_settings["spec"]
    )


def generate_fw(
    atoms, task_list, fw_settings, qadapter, update_settings=None, update_in_spec=True
):
    """Generates a FireWork

    Parameters
    ----------
    atoms: ase.atoms.Atoms, dict
        ASE Atoms object to preform the calculation on
    task_list: list of TaskSpecs
        Definitions for the tasks to be run
    fw_settings: dict
        FireWork settings for the step
    qadapter: dict
        The queueadapter for the step
    update_settings: dict
        update calculator settings
    update_in_spec: bool
        If True move the current out_spec to be in_spec

    Returns
    -------
    Firework
        A firework for the task
    """
    fw_settings = update_fw_settings(
        fw_settings, fw_settings["fw_name"], qadapter, update_in_spec=update_in_spec
    )
    fw_settings[
        "fw_name"
    ] += f"_{atoms.symbols.get_chemical_formula()}_{hash_atoms_and_calc(atoms)[0]}"

    if not update_settings:
        update_settings = {}

    if "in_spec_atoms" in fw_settings:
        at = fw_settings["in_spec_atoms"]
    else:
        at = atoms.copy()

    if "in_spec_calc" in fw_settings:
        cl = fw_settings["in_spec_calc"]
    else:
        cl = atoms.calc

    return generate_firework(
        task_list, at, cl, fw_settings, update_calc_settings=update_settings
    )


def generate_kgrid_fw(workflow, atoms, fw_settings):
    """Generate a k-grid optimization Firework

    Parameters
    ----------
    workflow: Settings
        workflow settings where the task is defined
    atoms: ase.atoms.Atoms, dict
        ASE Atoms object to preform the calculation on
    fw_settings: dict
        Firework settings for the step

    Returns
    -------
    Firework
        Firework for the k-grid optimization
    """
    # Get queue adapter settings
    fw_settings["fw_name"] = "kgrid_opt"
    fw_settings["out_spec_k_den"] = "kgrid"

    if "kgrid_qadapter" in workflow:
        qadapter = workflow["kgrid_qadapter"]
    else:
        qadapter = None

    func_kwargs = {
        "workdir": f"{workflow.general.workdir_cluster}/{fw_settings['fw_name']}/",
        "trajectory_file": filenames.trajectory,
        "dfunc_min": workflow.general.get("kgrid_dfunc_min", 1e-12),
    }

    if qadapter and "walltime" in qadapter:
        func_kwargs["walltime"] = str2time(qadapter["walltime"])
    else:
        func_kwargs["walltime"] = 1800

    task_spec = gen_kgrid_task_spec(func_kwargs)
    return generate_fw(atoms, task_spec, fw_settings, qadapter)


def generate_relax_fw(workflow, atoms, fw_settings, basisset_type):
    """Generates a Firework for the relaxation step

    Parameters
    ----------
    workflow: Settings
        workflow settings where the task is defined
    atoms: ase.atoms.Atoms, dict
        ASE Atoms object to preform the calculation on
    fw_settings: dict
        Firework settings for the step
    basisset_type: str
        Basis Set parameters to use for the calculation

    Returns
    -------
    Firework
        Firework for the relaxation step
    """
    if f"{basisset_type}_rel_qadapter" in workflow:
        qadapter = workflow["light_rel_qadapter"]
    else:
        qadapter = None

    abreviated_basis = [bt[0] for bt in basisset_type.split("_")]
    fw_settings["fw_name"] = f"{'_'.join(abreviated_basis)}_relax"

    func_kwargs = {
        "workdir": f"{workflow.general.workdir_cluster}/{fw_settings['fw_name']}/"
    }
    fw_out_kwargs = {"calc_number": 0}

    task_spec = gen_aims_task_spec(func_kwargs, fw_out_kwargs)

    if "relaxation" in workflow:
        method = workflow.relaxation.get("method", "trm")
        force_crit = workflow.relaxation.get("conv_crit", 1e-3)
        relax_unit_cell = workflow.relaxation.get("relax_unit_cell", "full")
    else:
        method = "trm"
        force_crit = 1e-3
        relax_unit_cell = "full"

    update_settings = {
        "relax_geometry": f"{method} {force_crit}",
        "basisset_type": basisset_type,
        "relax_unit_cell": relax_unit_cell,
    }

    return generate_fw(atoms, task_spec, fw_settings, qadapter, update_settings, True)


def generate_phonon_fw(workflow, atoms, fw_settings, typ):
    """Generates a Firework for the phonon initialization

    Parameters
    ----------
    aworkflow: Settings
        workflow settings where the task is defined
    atoms: ase.atoms.Atoms or dict
        ASE Atoms object to preform the calculation on
    fw_settings: dict
        Firework settings for the step
    typ: str
        either phonopy or phono3py

    Returns
    -------
    Firework
        Firework for the relaxation step
    """

    if f"{typ}_qadapter" in workflow:
        qadapter = workflow["phonopy_qadapter"]
    else:
        qadapter = {}

    update_settings = {}
    if "basisset_type" in workflow[typ]:
        update_settings["basisset_type"] = workflow[typ].pop("basisset_type")

    if "socket_io_port" in workflow[typ]:
        update_settings["use_pimd_wrapper"] = workflow[typ].pop("socket_io_port")
    elif "use_pimd_wrapper" in workflow[typ]:
        update_settings["use_pimd_wrapper"] = workflow[typ].pop("use_pimd_wrapper")

    if (
        workflow[typ].get("serial", True)
        and "spec" in fw_settings
        and "prev_dos_fp" in fw_settings["spec"]
    ):
        phonon, _, scs = preprocess(atoms, workflow[typ]["supercell_matrix"])
        qadapter["walltime"] = time2str(str2time(qadapter["walltime"]) * len(scs))
        if len(atoms) * np.linalg.det(phonon.get_supercell_matrix()) > 200:
            update_settings["use_local_index"] = True
            update_settings["load_balancing"] = True

    if "walltime" in qadapter:
        workflow[typ]["walltime"] = str2time(qadapter["walltime"])
    else:
        workflow[typ]["walltime"] = 1800

    fw_settings["fw_name"] = typ
    natoms = len(atoms) * np.linalg.det(
        get_3x3_matrix(workflow[typ]["supercell_matrix"])
    )
    natoms = int(round(natoms))
    workflow[typ][
        "workdir"
    ] = f"{workflow.general.workdir_cluster}/sc_natoms_{natoms}/{typ}/"
    if typ == "phonopy":
        func_kwargs = {"ph_settings": workflow[typ].copy()}
    elif typ == "phono3py":
        func_kwargs = {"ph3_settings": workflow[typ].copy()}

    task_spec = gen_phonon_task_spec(func_kwargs, fw_settings)

    return generate_fw(atoms, task_spec, fw_settings, qadapter, update_settings)


def generate_phonon_postprocess_fw(workflow, atoms, fw_settings, typ):
    """Generates a Firework for the phonon analysis

    Parameters
    ----------
    atoms: ase.atoms.Atoms, dict
        ASE Atoms object to preform the calculation on
    wd: str
        Workdirectory
    fw_settings: dict
        Firework settings for the step
    ph_settings: dict
        kwargs for the phonon analysis
    wd_init: str
        workdir for the initial phonon force calculations

    Returns
    -------
    Firework
        Firework for the phonon analysis
    """
    if typ == "phonopy":
        fw_settings["mod_spec_add"] = "ph"
        fw_settings["fw_name"] = "phonopy_analysis"
    else:
        fw_settings["fw_name"] = "phono3py_analysis"
        fw_settings["mod_spec_add"] = "ph3"
    fw_settings["mod_spec_add"] += "_forces"

    func_kwargs = workflow[typ].copy()
    if "workdir" in func_kwargs:
        func_kwargs.pop("workdir")

    natoms = len(atoms) * np.linalg.det(
        get_3x3_matrix(workflow[typ]["supercell_matrix"])
    )
    natoms = int(round(natoms))

    func_kwargs[
        "analysis_workdir"
    ] = f"{workflow.general.workdir_local}/sc_natoms_{natoms}/{fw_settings['fw_name']}/"
    func_kwargs["init_workdir"] = f"{workflow.general.workdir_cluster}/{typ}/"

    task_spec = gen_phonon_analysis_task_spec(
        "vibes." + fw_settings["fw_name"][:-9] + ".postprocess.postprocess",
        func_kwargs,
        fw_settings["mod_spec_add"][:-7] + "_metadata",
        fw_settings["mod_spec_add"],
        fw_settings["mod_spec_add"][:-7] + "_times",
        False,
    )
    fw_settings[
        "fw_name"
    ] += f"_{atoms.symbols.get_chemical_formula()}_{hash_atoms_and_calc(atoms)[0]}"
    return generate_firework(task_spec, None, None, fw_settings=fw_settings.copy())


def generate_phonon_fw_in_wf(
    atoms, wd, fw_settings, qadapter, ph_settings, update_in_spec=True
):
    """Generates a Firework for the phonon initialization

    Parameters
    ----------
    atoms: ase.atoms.Atoms, dict
        ASE Atoms object to preform the calculation on
    wd: str
        Workdirectory
    fw_settings: dict
        Firework settings for the step
    qadapter: dict
        The queueadapter for the step
    ph_settings: dict
        kwargs for the phonons
    update_settings: dict
        calculator update settings

    Returns
    -------
    Firework
        Firework for the phonon initialization
    """
    if (
        "serial" in ph_settings
        and ph_settings["serial"]
        and "spec" in fw_settings
        and "prev_dos_fp" in fw_settings["spec"]
    ):
        _, _, scs = preprocess(atoms, ph_settings["supercell_matrix"])
        qadapter["walltime"] = time2str(str2time(qadapter["walltime"]) * len(scs))

    if qadapter and "walltime" in qadapter:
        ph_settings["walltime"] = str2time(qadapter["walltime"])
    else:
        ph_settings["walltime"] = 1800

    update_settings = {}
    if "basisset_type" in ph_settings:
        update_settings["basisset_type"] = ph_settings.pop("basisset_type")
    if "socket_io_port" in ph_settings:
        update_settings["use_pimd_wrapper"] = ph_settings.pop("socket_io_port")
    elif "use_pimd_wrapper" in ph_settings:
        update_settings["use_pimd_wrapper"] = ph_settings.pop("use_pimd_wrapper")

    typ = ph_settings.pop("type")
    fw_settings["fw_name"] = typ
    ph_settings["workdir"] = wd + "/" + typ + "/"
    if typ == "phonopy":
        func_kwargs = {"ph_settings": ph_settings.copy()}
    else:
        func_kwargs = {"ph3_settings": ph_settings.copy()}
    task_spec = gen_phonon_task_spec(func_kwargs, fw_settings)
    return generate_fw(
        atoms, task_spec, fw_settings, qadapter, update_settings, update_in_spec
    )


def generate_phonon_postprocess_fw_in_wf(
    atoms, wd, fw_settings, ph_settings, wd_init=None
):
    """Generates a Firework for the phonon analysis

    Parameters
    ----------
    atoms: ase.atoms.Atoms or dict
        ASE Atoms object to preform the calculation on
    wd: str
        Workdirectory
    fw_settings: dict
        Firework settings for the step
    ph_settings: dict
        kwargs for the phonon analysis
    wd_init: str
        workdir for the initial phonon force calculations

    Returns
    -------
    Firework
        Firework for the phonon analysis
    """
    if ph_settings.pop("type") == "phonopy":
        fw_settings["mod_spec_add"] = "ph"
        fw_settings["fw_name"] = "phonopy_analysis"
    else:
        fw_settings["fw_name"] = "phono3py_analysis"
        fw_settings["mod_spec_add"] = "ph3"
    fw_settings["mod_spec_add"] += "_forces"

    func_kwargs = ph_settings.copy()
    if "workdir" in func_kwargs:
        func_kwargs.pop("workdir")
    func_kwargs["analysis_workdir"] = wd + "/" + fw_settings["fw_name"] + "/"
    func_kwargs["init_workdir"] = wd_init
    task_spec = gen_phonon_analysis_task_spec(
        "vibes." + fw_settings["fw_name"][:-9] + ".postprocess.postprocess",
        func_kwargs,
        fw_settings["mod_spec_add"][:-7] + "_metadata",
        fw_settings["mod_spec_add"],
        fw_settings["mod_spec_add"][:-7] + "_times",
        False,
    )
    fw_settings[
        "fw_name"
    ] += f"_{atoms.symbols.get_chemical_formula()}_{hash_atoms_and_calc(atoms)[0]}"
    return generate_firework(task_spec, None, None, fw_settings=fw_settings.copy())


def generate_stat_samp_fw(workflow, atoms, fw_settings):
    """
    Generates a Firework for the statistical sampling initialization
    Args:
        workflow (settings.Settings): workflow settings object
        atoms (ase.Atoms or dict): ASE Atoms object to preform the calculation on
        fw_settings (dict): Firework settings for the step

    Returns:
        (Firework): Firework for the harmonic analysis initialization
    """
    fw_settings["fw_name"] = "stat_samp"

    if "statistical_sampling_qadapter" in workflow:
        qadapter = workflow["statistical_sampling_qadapter"]
    elif "phonopy_qadapter" in workflow:
        qadapter = workflow["phonopy_qadapter"]
    else:
        qadapter = None

    if qadapter and "walltime" in qadapter:
        workflow.statistical_sampling["walltime"] = str2time(qadapter["walltime"])
    else:
        workflow.statistical_sampling["walltime"] = 1800

    add_qadapter = False
    if "phonopy" in workflow:
        add_qadapter = workflow.phonopy.get("converge_phonons", False)

    workflow.statistical_sampling[
        "workdir"
    ] = f"{workflow.general.workdir_cluster}/statistical_sampling/"
    if "phonon_file" not in workflow.statistical_sampling:
        if "phonopy" not in workflow:
            raise IOError("phonon file must be given")

        if workflow.phonopy.get("converge_phonons", False):
            workflow.statistical_sampling[
                "phonon_file"
            ] = f"{workflow.general.workdir_local}/converged/trajectory.son"
        else:
            sc_mat = get_3x3_matrix(workflow.phonopy.supercell_matrix)
            sc_natoms = int(round(np.linalg.det(sc_mat) * len(atoms)))
            rel_dir = f"/sc_natoms_{sc_natoms}/phonopy_analysis/trajectory.son"
            workflow.statistical_sampling["phonon_file"] = (
                workflow.general.workdir_local + rel_dir
            )
    fw_settings.pop("in_spec_calc", None)
    fw_settings.pop("in_spec_atoms", None)
    fw_settings["from_db"] = False

    task_spec = gen_stat_samp_task_spec(
        workflow.statistical_sampling, fw_settings, add_qadapter
    )
    return generate_fw(atoms, task_spec, fw_settings, qadapter, None, False)


def generate_stat_samp_postprocess_fw(workflow, atoms, fw_settings):
    """
    Generates a Firework for the statistical sampling analysis

    Args:
        workflow (settings.Settings): workflow settings object
        atoms (ase.Atoms or dict): ASE Atoms object to preform the calculation on
        fw_settings (dict): Firework settings for the step

    Returns
    -------
    Firework
        Firework for the anharmonicity analysis
    """
    fw_settings["fw_name"] = "statistical_sampling_analysis"
    fw_settings["mod_spec_add"] = "stat_samp"
    fw_settings["mod_spec_add"] += "_forces"

    func_kwargs = workflow["statistical_sampling"].copy()
    if "workdir" in func_kwargs:
        func_kwargs.pop("workdir")

    func_kwargs[
        "analysis_workdir"
    ] = f"{workflow.general.workdir_local}/{fw_settings['fw_name']}/"

    task_spec = gen_stat_samp_analysis_task_spec(
        func_kwargs,
        fw_settings["mod_spec_add"][:-7] + "_metadata",
        fw_settings["mod_spec_add"],
        False,
    )
    fw_settings[
        "fw_name"
    ] += f"_{atoms.symbols.get_chemical_formula()}_{hash_atoms_and_calc(atoms)[0]}"
    return generate_firework(task_spec, None, None, fw_settings=fw_settings.copy())


def generate_aims_fw(workflow, atoms, fw_settings):
    """Generates a Firework for the relaxation step

    Parameters
    ----------
    workflow: Settings
        workflow settings where the task is defined
    atoms: ase.atoms.Atoms or dict
        ASE Atoms object to preform the calculation on
    fw_settings: dict
        Firework settings for the step

    Returns
    -------
    Firework
        Firework for the relaxation step
    """
    if f"aims_qadapter" in workflow:
        qadapter = workflow["aims_qadapter"]
    else:
        qadapter = None

    fw_settings["fw_name"] = f"aims"

    func_kwargs = {"workdir": f"{workflow.general.workdir_cluster}/aims_calculation/"}
    task_spec = gen_aims_task_spec(func_kwargs, {}, relax=False)

    return generate_fw(atoms, task_spec, fw_settings, qadapter, None, True)


def generate_gruniesen_fd_fw(
    workflow, atoms, trajectory_file, constraints, fw_settings
):
    """Generate a FireWork to calculate the Gruniesen Parameter with finite differences

    Parameters
    ----------
    workflow: Settings
        The Workflow Settings
    atoms: ase.atoms.Atoms
        The initial ASE Atoms object of the primitive cell
    trajecotry_file: str
        Path the the equilibrium phonon trajectory
    constraints: list of dict
        list of relevant constraint dictionaries for relaxations
    fw_settings: dict
        The FireWork Settings for the current job

    Returns
    -------
    Firework:
        The Gruniesen setup firework
    """
    chem_form = atoms.symbols.get_chemical_formula(empirical=True, mode="metal")
    atoms_hash = hash_atoms_and_calc(atoms)[0]
    fw_settings["fw_name"] = f"gruniesen_setup_{chem_form}_{atoms_hash}"

    task_spec = gen_gruniesen_task_spec(workflow, trajectory_file, constraints)
    return generate_firework(task_spec, None, None, fw_settings.copy())


def generate_md_fw(workflow, atoms, fw_settings, qadapter=None, workdir=None):
    """Generate a FireWork to run a Molecular Dynamics calculation

    Parameters
    ----------
    workflow: Settings
        The Workflow Settings
    atoms: ase.atoms.Atoms
        The initial ASE Atoms object of the primitive cell
    fw_settings: dict
        The FireWork Settings for the current job
    qadapter: dict
        The queueadapter for the step
    workdir: str
        The working directory for the calculation

    Returns
    -------
    Firework:
        The Molecular Dynamics setup firework
    """
    fw_settings = fw_settings.copy()
    fw_settings.pop("in_spec_atoms", None)
    fw_settings.pop("in_spec_calc", None)
    fw_settings["from_db"] = False

    if qadapter is None:
        qadapter = workflow.pop("md_qadapter", None)

    md_settings = workflow["md"].copy()

    if "phonon_file" not in workflow.md and "phonopy" in workflow:
        if workflow.phonopy.get("converge_phonons", False):
            md_settings[
                "phonon_file"
            ] = f"{workflow.general.workdir_cluster}/converged/trajectory.son"
        else:
            sc_mat = get_3x3_matrix(workflow.phonopy.supercell_matrix)
            sc_natoms = int(round(np.linalg.det(sc_mat) * len(atoms)))
            rel_dir = f"/sc_natoms_{sc_natoms}/phonopy/trajectory.son"
            md_settings["phonon_file"] = workflow.general.workdir_cluster + rel_dir

    temps = md_settings.pop("temperatures", None)
    if temps is None:
        temps = [md_settings.pop("temperature")]

    fireworks = []
    for temp in temps:
        fw_settings["fw_name"] = f"md_{temp}"
        md_set = md_settings.copy()
        if workdir is None:
            md_set[
                "workdir"
            ] = f"{workflow.general.workdir_cluster}/{fw_settings['fw_name']}/"
        else:
            md_set["workdir"] = workdir
        md_set["temperature"] = temp
        task_spec = gen_md_task_spec(md_set, fw_settings)
        fireworks.append(
            generate_fw(atoms, task_spec, fw_settings, qadapter, None, False)
        )
    return fireworks
