"""Functions used to wrap around HiLDe Phonopy/Phono3py functions"""
from pathlib import Path

from hilde.helpers.converters import input2dict
from hilde.phonon_db.ase_converters import dict2atoms
from hilde.phonopy import displacement_id_str
from hilde.phonopy.workflow import bootstrap
from hilde.settings import Settings, AttributeDict
from hilde.statistical_sampling.workflow import bootstrap as bootstrap_stat_samp
from hilde.fireworks.tasks.general_py_task import get_func
from hilde.trajectory import step2file, metadata2file


def setup_calc(settings, calc, use_pimd_wrapper, kwargs_boot):
    """
    Sets up a calculation
    Parameters:
        settings (Settings): The settings object for the calculation
        calc (ASE Calculator): Calculator used for the calculation
        use_pimd_wrapper (dict): Dictionary to wrapper ipi parameters for calc
        kwargs_boot (dict): kwargs for the bootstrapping

    Returns:
        settings (Settings): The updated settings object
        kwargs_boot (dict): The updated kwargs for the bootstrapping

    """
    if calc.name.lower() != "aims":
        kwargs_boot["calculator"] = calc
    else:
        settings["control"] = calc.parameters.copy()
        if "species_dir" in settings.control:
            sd = settings["control"].pop("species_dir")
            settings["basisset"] = AttributeDict({"type": sd.split("/")[-1]})
        if use_pimd_wrapper:
            settings["socketio"] = AttributeDict(
                {"port": settings.pop("use_pimd_wrapper")}
            )
        if "aims_command" in settings.control:
            del settings.control["aims_command"]
    return settings, kwargs_boot


def setup_phonon_outs(ph_settings, settings, prefix, atoms, calc):
    """
    Sets up the phonon outputs
    Parameters:
        ph_settings (dict): Settings object for the phonopy, phono3py object
        settings (Settings): General settings for the step
        prefix (str): key prefix for the task
        atoms (ASE Atoms): ASE Atoms object for the material
        calc (ASE Calculator): Calculator used for the calculation

    Returns:
        out (dict): All the necessary output/metadata for the task
    """
    settings, kwargs_boot = setup_calc(
        settings,
        calc,
        ("use_pimd_wrapper" in settings and settings["use_pimd_wrapper"]),
        dict(),
    )
    settings[f"{prefix}onopy"] = ph_settings.copy()
    if "serial" in settings[f"{prefix}onopy"]:
        del settings[f"{prefix}onopy"]["serial"]
    out = bootstrap(name=f"{prefix}onopy", settings=settings, **kwargs_boot)

    out["metadata"]["supercell"] = {"atoms": out["metadata"]["atoms"], "calculator": {}}
    out["metadata"]["primitive"] = input2dict(atoms)
    out["prefix"] = prefix
    out["settings"] = ph_settings.copy()
    return out


def bootstrap_phonon(
    atoms, calc, kpt_density=None, ph_settings=None, ph3_settings=None, fw_settings=None
):
    """
    Creates a Settings object and passes it to the bootstrap function
    Parameters:
        atoms (ASE Atoms Object): Atoms object of the primitive cell
        calc (ASE Calculator): Calculator for the force calculations
        kpt_density (float): k-point density for the MP-Grid
        ph_settings (dict): kwargs for phonopy setup
        ph3_settings (dict): kwargs for phono3py setup
        fw_settings (dict): FireWork specific settings

    Returns:
        (dict): The output of hilde.phonopy.workflow.bootstrap for phonopy and phono3py
    """
    settings = Settings(settings_file=None)
    settings.atoms = atoms
    if kpt_density:
        settings["control_kpt"] = AttributeDict({"density": kpt_density})

    outputs = []
    at = atoms.copy()
    at.set_calculator(None)
    if ph_settings:
        outputs.append(setup_phonon_outs(ph_settings, settings, "ph", at, calc))

    if ph3_settings:
        outputs.append(setup_phonon_outs(ph3_settings, settings, "ph3", at, calc))
    return outputs


def bootstrap_stat_sample(
    atoms, calc, kpt_density=None, stat_samp_settings=None, fw_settings=None
):
    """
    Initializes the statistical sampling task
    Parameters:
        atoms (ASE Atoms Object): Atoms object of the primitive cell
        calc (ASE Calculator): Calculator for the force calculations
        kpt_density (float): k-point density for the MP-Grid
        stat_samp_settings (dict): kwargs for statistical sampling setup
        fw_settings (dict): FireWork specific settings

    Returns:
        (dict): The output of hilde.statistical_sampling.workflow.bootstrap
    """
    settings = Settings(settings_file=None)
    settings.atoms = atoms
    if kpt_density:
        settings["control_kpt"] = AttributeDict({"density": kpt_density})

    outputs = []
    at = atoms.copy()
    at.set_calculator(None)
    if stat_samp_settings:
        settings, kwargs_boot = setup_calc(
            settings,
            calc,
            (
                "use_pimd_wrapper" in stat_samp_settings
                and stat_samp_settings["use_pimd_wrapper"]
            ),
            dict(),
        )
        settings["statistical_sampling"] = stat_samp_settings.copy()
        stat_samp_out = bootstrap_stat_samp(
            name="statistical_sampling", settings=settings, **kwargs_boot
        )
        stat_samp_out["prefix"] = "stat_samp"
        stat_samp_out["settings"] = stat_samp_settings.copy()
        outputs.append(stat_samp_out)
    else:
        raise IOError("Stat sampling requires a settings object")
    return outputs


def collect_to_trajectory(workdir, trajectory, calculated_atoms, metadata):
    """
    Collects forces to a single trajectory file
    Parameters:
        workdir (str): working directory for the task
        trajectory (str): file name for the trajectory file
        calculated_atoms (list of ASE Atoms): Results of the force calculations
        metadata (dict): metadata for the phonon calculations
    """
    traj = Path(workdir) / trajectory
    traj.parent.mkdir(exist_ok=True, parents=True)
    if "Phonopy" in metadata:
        for el in metadata["Phonopy"]["displacement_dataset"]["first_atoms"]:
            el["number"] = int(el["number"])

    if "Phono3py" in metadata:
        for el1 in metadata["Phono3py"]["displacement_dataset"]["first_atoms"]:
            el1["number"] = int(el1["number"])
            for el2 in el1["second_atoms"]:
                el2["number"] = int(el2["number"])

    metadata2file(metadata, str(traj))
    if isinstance(calculated_atoms[0], dict):
        temp_atoms = [dict2atoms(cell) for cell in calculated_atoms]
    else:
        temp_atoms = calculated_atoms.copy()
    calculated_atoms = sorted(
        temp_atoms,
        key=lambda x: x.info[displacement_id_str] if x else len(calculated_atoms) + 1,
    )
    for atoms in calculated_atoms:
        if atoms:
            step2file(atoms, atoms.calc, str(traj))


def phonon_postprocess(func_path, phonon_times, max_mem, **kwargs):
    """
    @brief      performs phonon postprocessing steps

    @param      func_path     The path to the postprocessing function
    @param      phonon_times  The time it took to calculate the phonon forces
    @param      max_mem       Maximum memory useage of the calculation
    @param      kwargs        The keyword arguments for the phonon calculations

    @return     { description_of_the_return_value }
    """
    func = get_func(func_path)
    return func(**kwargs)
