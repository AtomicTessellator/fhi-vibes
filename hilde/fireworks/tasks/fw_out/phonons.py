"""Generate FWActions after setting Phonon Calculations"""

from ase.symbols import Symbols
from fireworks import FWAction, Workflow

from hilde.fireworks.workflows.firework_generator import (
    generate_phonon_fw_in_wf,
    generate_phonon_postprocess_fw_in_wf,
    generate_firework,
    time2str,
)
from hilde.phonon_db.ase_converters import calc2dict, atoms2dict
from hilde.helpers.converters import dict2atoms
from hilde.helpers.k_grid import k2d
from hilde.helpers.watchdogs import str2time
from hilde.phonon_db.row import phonon_to_dict, phonon3_to_dict
from hilde.structure.convert import to_Atoms
from hilde.trajectory import reader
from hilde.fireworks.tasks.postprocess.phonons import get_converge_phonon_update


def post_init_mult_calcs(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """postprocessing for initializing parallel force calculaitons

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        atoms reference structure for the calculation
    calc: ase.calculators.calulator.Calculator
        The claculator of the claulation
    outputs: dict
        The outputs after setting up claculations
    func: str
        Function path to the main function
    func_fw_out: str
        Function path to the fw_out function
    func_kwargs: dict
        The kwargs for the main function
    func_fw_kwargs: dict
        The kwargs for the fw_out function
    fw_settings: dict
        The FireWorks settings

    Returns
    -------
    FWAction
        The action that will run all force calculations
    """
    if fw_settings is None:
        fw_settings = dict()
    update_spec = dict()
    detours = []

    for out in outputs:
        prefix = out["prefix"]
        fw_set = fw_settings.copy()
        func_set = out["settings"]

        if prefix + "_settings" in func_fw_kwargs:
            func_set = dict(func_set, **func_fw_kwargs[prefix + "_settings"])
        if "serial" not in func_set:
            func_set["serial"] = True

        update_spec[prefix + "_metadata"] = out["metadata"]

        if "k_pt_density" in out:
            update_spec["kgrid"] = out["k_pt_density"]

        if "spec" in fw_set:
            fw_set["spec"].update(update_spec)
        else:
            fw_set["spec"] = update_spec.copy()

        fw_set["mod_spec_add"] = prefix + "_forces"
        fw_set["metadata_spec"] = prefix + "_metadata"

        if "calculator" in out:
            calc_dict = calc2dict(out["calculator"])
        else:
            calc_dict = calc.copy()

        if (
            "prev_dos_fp" in fw_set
            and "walltime" in fw_set
            and "serial" in func_set
            and func_set["serial"]
        ):
            fw_set["walltime"] = time2str(
                str2time(fw_set["walltime"]) * len(out["atoms_to_calculate"])
            )

        detours = get_detours(
            out["atoms_to_calculate"],
            calc_dict,
            prefix,
            func_set,
            fw_set,
            update_spec,
            atoms,
            detours,
        )
    return FWAction(update_spec=update_spec, detours=detours)


def get_detours(
    atoms_to_calculate,
    calc_dict,
    prefix,
    calc_kwargs,
    fw_settings,
    update_spec,
    atoms=None,
    detours=None,
):
    """Add a set of detours for force calculations

    Parameters
    ----------
    atoms_to_calculate: list of ase.atoms.Atoms
        List of structures to calculate forces for
    calc_dict: dict
        Dictionary representation of the ase.calculators.calulator.Calculator
    prefix: str
        prefix to add to force calculations
    calc_kwargs: dict
        A set of kwargs for the Force calculations
    fw_settings: dict
        A dictionary describing all FireWorks settings
    update_spec: dict
        Parmeters to be added to the FireWorks spec
    atoms: ase.atoms.Atoms
        Initial ASE Atoms object representation of the structure
    detours: list of Fireworks
        Current list of force calculations to perform

    Returns
    -------
    list of Fireworks
        The updated detours object
    """
    if detours is None:
        detours = []
    fw_settings["time_spec_add"] = prefix + "_times"
    if "walltime" in calc_kwargs:
        if "spec" in fw_settings and "_queueadapter" in fw_settings["spec"]:
            if "walltime" in fw_settings["spec"]["_queueadapter"]:
                calc_kwargs["walltime"] = (
                    str2time(fw_settings["spec"]["_queueadapter"]["walltime"]) - 120
                )
            else:
                fw_settings["spec"]["_queueadapter"]["walltime"] = (
                    time2str(calc_kwargs["walltime"]) + 120
                )
        else:
            if "spec" not in fw_settings:
                fw_settings["spec"] = dict()
            fw_settings["spec"]["_queueadapter"] = {
                "walltime": time2str(calc_kwargs["walltime"] + 120)
            }

    if calc_kwargs["serial"]:
        update_spec[prefix + "_calculated_atoms"] = [
            atoms2dict(at) for at in atoms_to_calculate
        ]
        update_spec[prefix + "_calculator"] = calc_dict
        fw_settings["spec"].update(update_spec)
        fw_settings["calc_atoms_spec"] = prefix + "_calculated_atoms"
        fw_settings["calc_spec"] = prefix + "_calculator"
        return add_socket_calc_to_detours(
            detours, atoms, calc_kwargs, fw_settings, prefix
        )
    return add_single_calc_to_detours(
        detours, calc_kwargs, atoms, atoms_to_calculate, calc_dict, fw_settings, prefix
    )


def add_socket_calc_to_detours(detours, atoms, func_kwargs, fw_settings, prefix):
    """Generates a Firework to run a socket calculator and adds it to the detours

    Parameters
    ----------
    detours: list of Fireworks
        Current list of detours
    atoms: ase.atoms.Atoms
        Initial ASE Atoms object representation of the structure
    func_kwargs: dict
        kwargs needed to do the socket I/O calculation
    fw_settings: dict
        FireWorks settings
    prefix: str
        ph for phonopy and ph3 for phono3py calculations

    Returns
    -------
    list of Fireworks
        an updated detours list
    """
    calc_kwargs = {}
    calc_keys = ["trajectory", "workdir", "backup_folder", "walltime"]
    for key in calc_keys:
        if key in func_kwargs:
            calc_kwargs[key] = func_kwargs[key]
    fw_set = fw_settings.copy()
    fw_set["fw_name"] = (
        prefix + f"_serial_forces_{Symbols(atoms['numbers']).get_chemical_formula()}"
    )
    fw = generate_firework(
        func="hilde.fireworks.tasks.calculate_wrapper.wrap_calc_socket",
        func_fw_out="hilde.fireworks.tasks.fw_out.calculate.socket_calc_check",
        func_kwargs=calc_kwargs,
        atoms_calc_from_spec=False,
        inputs=[
            prefix + "_calculated_atoms",
            prefix + "_calculator",
            prefix + "_metadata",
        ],
        fw_settings=fw_set.copy(),
    )
    detours.append(fw)
    return detours


def add_single_calc_to_detours(
    detours, func_fw_kwargs, atoms, atoms_list, calc_dict, fw_settings, prefix
):
    """Adds a group of Fireworks to run as single calculations

    Parameters
    ----------
    detours: list of Fireworks
        Current list of detours
    func_kwargs: dict
        kwargs needed to do the socket I/O calculation
    atoms: dict
        Dictionary representing the ASE Atoms object of theprimitive cell
    atoms_list: list of Atoms
        List of supercells to perform force calculations on
    calc_dict: dict
        Dictionary representing the ASE Calculator for the force calculations
    fw_settings: dict
        FireWorks settings
    prefix: str
        ph for phonopy and ph3 for phono3py calculations

    Returns
    -------
    list of Fireworks
        an updated detours list
    """
    for i, sc in enumerate(atoms_list):
        if not sc:
            continue
        fw_settings = fw_settings.copy()
        fw_settings["from_db"] = False
        fw_settings.pop("kpoint_density_spec", None)
        sc.info["displacement_id"] = i
        sc_dict = atoms2dict(sc)
        for key, val in calc_dict.items():
            sc_dict[key] = val
        calc_kwargs = {"workdir": func_fw_kwargs["workdir"] + f"/{i:05d}"}
        fw_settings["fw_name"] = (
            prefix + f"forces_{Symbols(sc_dict['numbers']).get_chemical_formula()}_{i}"
        )
        detours.append(
            generate_firework(
                func="hilde.fireworks.tasks.calculate_wrapper.wrap_calculate",
                func_fw_out="hilde.fireworks.tasks.fw_out.calculate.mod_spec_add",
                func_kwargs=calc_kwargs,
                atoms=sc_dict,
                calc=calc_dict,
                atoms_calc_from_spec=False,
                fw_settings=fw_settings,
            )
        )
    return detours


def add_phonon_to_spec(func, func_fw_out, *args, fw_settings=None, **kwargs):
    """Add the phonon_dict to the spec

    Parameters
    ----------
    func: str
        Path to the phonon analysis function
    func_fw_out: str
        Path to this function
    fw_settings: dict
        Dictionary for the FireWorks specific systems
    kwargs: dict
        Dictionary of keyword arguments that must have the following objects

        ouputs: phonopy.Phonopy
            The Phonopy object from post-processing

    Returns
    -------
    FWAction
        FWAction that adds the phonon_dict to the spec
    """
    traj = f"{kwargs['workdir']}/{kwargs['trajectory']}"
    _, metadata = reader(traj, True)
    calc_dict = metadata["calculator"]
    calc_dict["calculator"] = calc_dict["calculator"].lower()
    calc_dict["calculator"] == "aims"
    if calc_dict["calculator"] == "aims":
        k_pt_density = k2d(
            dict2atoms(metadata["supercell"]["atoms"]),
            calc_dict["calculator_parameters"]["k_grid"],
        )
    else:
        k_pt_density = None
    qadapter = dict()
    if fw_settings and "spec" in fw_settings:
        qadapter = fw_settings["spec"].get("_queueadapter", None)
    if "phono3py" in args[0]:
        update_spec = {
            "ph3_dict": phonon3_to_dict(kwargs["outputs"]),
            "ph3_calculator": calc_dict,
            "ph3_supercell": atoms2dict(to_Atoms(kwargs["outputs"].get_primitive())),
        }
    else:
        update_spec = {
            "ph_dict": phonon_to_dict(kwargs["outputs"]),
            "ph_calculator": calc_dict,
            "ph_supercell": atoms2dict(to_Atoms(kwargs["outputs"].get_primitive())),
            "_queueadapter": qadapter,
        }
    update_spec["kgrid"] = k_pt_density
    return FWAction(update_spec=update_spec)


def converge_phonons(func, func_fw_out, *args, fw_settings=None, **kwargs):
    """Check phonon convergence and set up future calculations after a phonon calculation

    Parameters
    ----------
    func: str
        Path to the phonon analysis function
    func_fw_out: str
        Path to this function
    args: list
        list arguments passed to the phonon analysis
    fw_settings: dict
        Dictionary for the FireWorks specific systems
    kwargs: dict
        Dictionary of keyword arguments with the following keys

        outputs: phonopy.Phonopy
            The Phonopy object from post-processing
        serial: bool
            If True use a serial calculation
        init_workdir: str
            Path to the base phonon force calculations
        trajectory: str
            trajectory file name

    Returns
    -------
    FWAction
        Increases the supercell size or adds the phonon_dict to the spec
    """
    if not fw_settings:
        fw_settings = dict()
        fw_settings["spec"] = dict()
    elif "spec" not in fw_settings:
        fw_settings["spec"] = dict()

    fw_settings["from_db"] = False
    fw_settings.pop("in_spec_calc", None)
    fw_settings.pop("in_spec_atoms", None)

    fw_settings["spec"]["kgrid"] = args[-1]

    ph = kwargs.pop("outputs")
    workdir = kwargs.pop("workdir")
    trajectory = kwargs.pop("trajectory")
    conv, update_job = get_converge_phonon_update(
        workdir, trajectory, args[1], ph, **kwargs
    )

    if conv:
        qadapter = None
        qadapter = fw_settings["spec"].get("_queueadapter", None)
        update_spec = dict(_queueadapter=qadapter, **update_job)
        update_spec["kgrid"] = args[-1]
        return FWAction(update_spec=update_spec)

    fw_settings["in_spec_calc"] = "ph_calculator"
    update_spec = {
        "ph_calculator": update_job["ph_calculator"],
        "prev_dos_fp": update_job["prev_dos_fp"],
    }

    fw_settings["spec"].update(update_spec)

    update_spec["ph_supercell"] = update_job["ph_supercell"]
    func_kwargs = {
        "type": "phonopy",
        "displacement": update_job["displacement"],
        "supercell_matrix": update_job["supercell_matrix"],
        "serial": kwargs["serial"],
        "converge_phonons": True,
        "prev_wd": update_job["prev_wd"],
    }
    kwargs.update(func_kwargs)

    fw_settings["spec"]["prev_dos_fp"] = update_job["prev_dos_fp"]

    if "_queueadapter" in fw_settings["spec"]:
        func_kwargs.pop("walltime", None)
        fw_settings["spec"]["_queueadapter"].pop("queue", None)
        fw_settings["spec"]["_queueadapter"]["walltime"] = time2str(
            update_job["expected_walltime"]
        )
        # fw_settings["spec"]["_queueadapter"]["expected_mem"] = update_job[
        #     "expected_mem"
        # ]
        fw_settings["spec"]["_queueadapter"]["ntasks"] = update_job["ntasks"]
        qadapter = fw_settings["spec"]["_queueadapter"]
    else:
        qadapter = None

    pc = to_Atoms(ph.get_primitive())
    init_fw = generate_phonon_fw_in_wf(
        pc,
        update_job["init_workdir"],
        fw_settings,
        qadapter,
        func_kwargs,
        update_in_spec=False,
    )

    kwargs["prev_dos_fp"] = update_job["prev_dos_fp"]
    kwargs["trajectory"] = trajectory.split("/")[-1]
    kwargs["sc_matrix_original"] = update_job["sc_matrix_original"]
    analysis_fw = generate_phonon_postprocess_fw_in_wf(
        pc,
        update_job["analysis_wd"],
        fw_settings,
        kwargs,
        wd_init=update_job["init_workdir"],
    )

    detours = [init_fw, analysis_fw]
    wf = Workflow(detours, {init_fw: [analysis_fw]})
    return FWAction(detours=wf, update_spec=update_spec)
