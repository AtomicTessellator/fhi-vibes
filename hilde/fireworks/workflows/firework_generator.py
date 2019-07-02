"""Generates Task Specific FireWorks"""
import numpy as np

from fireworks import Firework

from hilde.fireworks.workflows.task_generator import (
    generate_task,
    generate_update_calc_task,
    generate_mod_calc_task,
)
from hilde.fireworks.tasks.task_spec import TaskSpec
from hilde.fireworks.tasks.utility_tasks import update_calc
from hilde.fireworks.workflows.task_spec_generator import (
    gen_phonon_task_spec,
    gen_stat_samp_task_spec,
    gen_phonon_analysis_task_spec,
    gen_aims_task_spec,
    gen_kgrid_task_spec,
)
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.helpers.k_grid import update_k_grid, update_k_grid_calc_dict
from hilde.phonon_db.ase_converters import atoms2dict, calc2dict
from hilde.phonopy.wrapper import preprocess


def get_time(time_str):
    """Converts a time step to the number of seconds that time stamp represents

    Parameters
    ----------
    time_str: str
        A string representing a specified time

    Returns
    -------
    n_secs: int
        A time presented as a number of seconds
    """
    time_set = time_str.split(":")
    n_secs = 0
    for ii, tt in enumerate(time_set):
        n_secs += int(round(float(tt))) * 60 ** (len(time_set) - 1 - ii)
    return int(n_secs)


def to_time_str(n_sec):
    """Converts a number of seconds into a time string

    Parameters
    ----------
    n_secs: int
        A time presented as a number of seconds

    Returns
    -------
    time_str: str
        A string representing a specified time
    """
    secs = int(n_sec % 60)
    mins = int(n_sec / 60) % 60
    hrs = int(n_sec / 3600)
    return f"{hrs}:{mins}:{secs}"


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
    calc=None,
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
    """A function that takes in a set of inputs and returns a Firework to perform that operation

    Parameters
    ----------
    task_spec_list: list of TaskSpecs
        list of task specifications to perform
    atoms: ase.atoms.Atoms, dictionary or str
        If not atoms_calc_from_spec then this must be an ASE Atoms object or a dictionary describing it If atoms_calc_from_spec then this must be a key str to retrieve the Atoms Object from the MongoDB launchpad
    calc: ase.calculators.calulator.Calculator, dictionary or str
        If not atoms_calc_from_spec then this must be an ase.calculators.calulator.Calculator or a dictionary describing it If atoms_calc_from_spec then this must be a key str to retrieve the Calculator from the MongoDB launchpad
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
        A Firework that will perform the desired operation on a set of atoms, and process the outputs for Fireworks

    Raises
    ------
    IOError
        If conflicting task_spec definitions provided, or none are provided
    """
    fw_settings = fw_settings.copy()
    if "spec" not in fw_settings:
        fw_settings["spec"] = dict()
    if update_calc_settings is None:
        update_calc_settings = dict()

    if func:
        if task_spec_list:
            raise IOError(
                "You have defined a task_spec_list and arguments to generate one, please only specify one of these"
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
            "You have not defined a task_spec_list or arguments to generate one, please specify one of these"
        )
    if isinstance(task_spec_list, TaskSpec):
        task_spec_list = [task_spec_list]

    if fw_settings and "from_db" in fw_settings:
        atoms_calc_from_spec = fw_settings["from_db"]

    if "fw_name" not in fw_settings:
        fw_settings["fw_base_name"] = ""
    elif "fw_base_name" not in fw_settings:
        fw_settings["fw_base_name"] = fw_settings["fw_name"]

    setup_tasks = []
    if atoms:
        if not atoms_calc_from_spec:
            # Preform calc updates here
            at = atoms2dict(atoms)
            if not isinstance(calc, str):
                if "k_grid_density" in update_calc_settings:
                    if not isinstance(calc, dict):
                        update_k_grid(
                            atoms, calc, update_calc_settings["k_grid_density"]
                        )
                    else:
                        recipcell = np.linalg.pinv(at["cell"]).transpose()
                        calc = update_k_grid_calc_dict(
                            calc,
                            recipcell,
                            at["pbc"],
                            update_calc_settings["k_grid_density"],
                        )

                cl = calc2dict(calc)

                for key, val in update_calc_settings.items():
                    if key != "k_grid_density":
                        cl = update_calc(cl, key, val)
                for key, val in cl.items():
                    at[key] = val
            else:
                cl = calc
                setup_tasks.append(
                    generate_update_calc_task(calc, update_calc_settings)
                )
        else:
            # Add tasks to update calculator parameters
            at = atoms
            cl = calc
            if update_calc_settings.keys():
                setup_tasks.append(
                    generate_update_calc_task(calc, update_calc_settings)
                )

        if "kpoint_density_spec" in fw_settings:
            setup_tasks.append(
                generate_mod_calc_task(
                    at, cl, "calculator", fw_settings["kpoint_density_spec"]
                )
            )
            cl = "calculator"
    else:
        at = None
        cl = None
    job_tasks = []
    for task_spec in task_spec_list:
        job_tasks.append(generate_task(task_spec, fw_settings, at, cl))
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

    # Get convergence criteria
    if "kgrid_dfunc_min" in workflow.general:
        dfunc_min = workflow.general.kgrid_dfunc_min
    else:
        dfunc_min = 1e-12

    func_kwargs = {
        "workdir": f"{workflow.general.workdir_cluster}/{fw_settings['fw_name']}/",
        "trajectory": "trajectory.son",
        "dfunc_min": dfunc_min,
    }

    if qadapter and "walltime" in qadapter:
        func_kwargs["walltime"] = get_time(qadapter["walltime"])
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

    fw_settings["fw_name"] = f"{basisset_type}_relax"

    func_kwargs = {
        "workdir": f"{workflow.general.workdir_cluster}/{fw_settings['fw_name']}/"
    }
    fw_out_kwargs = {"relax_step": 0}

    task_spec = gen_aims_task_spec(func_kwargs, fw_out_kwargs)

    if "relaxation" in workflow:
        method = workflow.relaxation.get("method", "trm")
        force_crit = workflow.relaxation.get("conv_crit", 1e-3)
    else:
        method = "trm"
        force_crit = 1e-3

    update_settings = {
        "relax_geometry": f"{method} {force_crit}",
        "relax_unit_cell": "full",
        "basisset_type": basisset_type,
        "scaled": True,
        "use_sym": True,
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
        qadapter = dict()

    if (
        workflow[typ].get("serial", True)
        and "spec" in fw_settings
        and "prev_dos_fp" in fw_settings["spec"]
    ):
        _, _, scs = preprocess(atoms, workflow[typ]["supercell_matrix"])
        qadapter["walltime"] = to_time_str(get_time(qadapter["walltime"]) * len(scs))

    if "walltime" in qadapter:
        workflow[typ]["walltime"] = get_time(qadapter["walltime"])
    else:
        workflow[typ]["walltime"] = 1800

    update_settings = {}
    if "basisset_type" in workflow[typ]:
        update_settings["basisset_type"] = workflow[typ].pop("basisset_type")

    if "socket_io_port" in workflow[typ]:
        update_settings["use_pimd_wrapper"] = workflow[typ].pop("socket_io_port")
    elif "use_pimd_wrapper" in workflow[typ]:
        update_settings["use_pimd_wrapper"] = workflow[typ].pop("use_pimd_wrapper")

    fw_settings["fw_name"] = typ
    workflow[typ]["workdir"] = f"{workflow.general.workdir_cluster}/{typ}/"
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

    func_kwargs[
        "analysis_workdir"
    ] = f"{workflow.general.workdir_local}/{fw_settings['fw_name']}/"
    func_kwargs["init_wd"] = f"{workflow.general.workdir_cluster}/{typ}/"

    task_spec = gen_phonon_analysis_task_spec(
        "hilde." + fw_settings["fw_name"][:-9] + ".postprocess.postprocess",
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
        qadapter["walltime"] = to_time_str(get_time(qadapter["walltime"]) * len(scs))

    if qadapter and "walltime" in qadapter:
        ph_settings["walltime"] = get_time(qadapter["walltime"])
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
    func_kwargs["init_wd"] = wd_init
    task_spec = gen_phonon_analysis_task_spec(
        "hilde." + fw_settings["fw_name"][:-9] + ".postprocess.postprocess",
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


# def generate_stat_samp_fw(workflow, atoms, fw_settings):
#     """Generates a Firework for the phonon initialization

#     Parameters
#     ----------
#     atoms: ase.atoms.Atoms, dict
#         ASE Atoms object to preform the calculation on
#     wd: str
#         Workdirectory
#     fw_settings: dict
#         Firework settings for the step
#     qadapter: dict
#         The queueadapter for the step
#     stat_samp_settings: dict
#         kwargs for the harmonic analysis

#     Returns
#     -------
#     Firework
#         Firework for the harmonic analysis initialization
#     """
#     if "statistical_sampling_qadapter" in workflow:
#         qadapter = workflow["statistical_sampling_qadapter"]
#     elif "phonopy_qadapter" in workflow:
#         qadapter = workflow["phonopy_qadapter"]
#     else:
#         qadapter = None

#     if qadapter and "walltime" in qadapter:
#         workflow.statistical_sampling["walltime"] = get_time(qadapter["walltime"])
#     else:
#         workflow.statistical_sampling["walltime"] = 1800

#     workflow.statistical_sampling[
#         "workdir"
#     ] = f"{workflow.general.workdir_cluster}/statistical_sampling/"
#     if workflow.phonopy.get("converge_phonons", False):
#         workflow.statistical_sampling[
#             "phonon_file"
#         ] = f"{workflow.general.workdir_local}/converged/trajectory.son"
#     else:
#         workflow.statistical_sampling[
#             "phonon_file"
#         ] = f"{workflow.general.workdir_local}/phonopy_analysis/trajectory.son"

#     fw_settings["fw_name"] = "statistical_sampling"
#     fw_settings["time_spec_add"] = "stat_samp_times"

#     task_spec = gen_stat_samp_task_spec(workflow.statistical_sampling, fw_settings)

#     return generate_fw(atoms, task_spec, fw_settings, qadapter, None, False)


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

    fw_settings["fw_name"] = f"aims_calculation"

    func_kwargs = {"workdir": f"{workflow.general.workdir_cluster}/aims_calculation/"}
    task_spec = gen_aims_task_spec(func_kwargs, dict(), relax=False)

    return generate_fw(atoms, task_spec, fw_settings, qadapter, None, True)
