"""Functions used to wrap around HiLDe Phonopy/Phono3py functions"""
from pathlib import Path

import numpy as np
from ase.constraints import (
    FixCartesianParametricRelations,
    FixScaledParametricRelations,
    dict2constraint,
)

from vibes.aims.context import AimsContext
from vibes.aims.setup import setup_aims
from vibes.fireworks.tasks.general_py_task import get_func
from vibes.fireworks.workflows.workflow_generator import generate_workflow
from vibes.helpers.converters import dict2atoms, input2dict
from vibes.helpers.numerics import get_3x3_matrix
from vibes.helpers.supercell import make_supercell
from vibes.phonopy import displacement_id_str
from vibes.phonopy.context import PhonopyContext
from vibes.phonopy.postprocess import postprocess
from vibes.phonopy.workflow import bootstrap
from vibes.settings import AttributeDict, Settings, TaskSettings
from vibes.structure.convert import to_Atoms
from vibes.trajectory import metadata2file, reader, step2file


def setup_calc(settings, calc, use_pimd_wrapper, kwargs_boot):
    """Sets up a calculation

    Parameters
    ----------
    settings: Settings
        The settings object for the calculation
    calc: ase.calculators.calulator.Calculator
        Calculator used for the calculation
    use_pimd_wrapper: dict
        Dictionary to wrapper ipi parameters for calc
    kwargs_boot: dict
        kwargs for the bootstrapping

    Returns
    -------
    settings: Settings
        The updated settings object
    kwargs_boot: dict
        The updated kwargs for the bootstrapping

    """
    if calc.name.lower() != "aims":
        kwargs_boot["calculator"] = calc
    else:
        settings["control"] = calc.parameters.copy()
        if "species_dir" in settings.control:
            sd = settings["control"].pop("species_dir")
            settings["basissets"] = AttributeDict({"default": sd.split("/")[-1]})
        if use_pimd_wrapper:
            settings["socketio"] = AttributeDict(
                {"port": settings.pop("use_pimd_wrapper")}
            )
        settings.control.pop("aims_command", None)
        if "control_kpt" in settings:
            settings.control.pop("kgrid", None)
    return settings, kwargs_boot


def setup_phonon_outputs(ph_settings, settings, prefix, atoms, calc):
    """Sets up the phonon outputs

    Parameters
    ----------
    ph_settings: dict
        Settings object for the phonopy, phono3py object
    settings: Settings
        General settings for the step
    prefix: str
        key prefix for the task
    atoms: ase.atoms.Atoms
        ASE Atoms object for the material
    calc: ase.calculators.calulator.Calculator
        Calculator used for the calculation

    Returns
    -------
    outputs: dict
        All the necessary output/metadata for the task
    """
    settings, kwargs_boot = setup_calc(
        settings,
        calc,
        ("use_pimd_wrapper" in settings and settings["use_pimd_wrapper"]),
        {},
    )
    settings[f"{prefix}onopy"] = ph_settings.copy()
    if "serial" in settings[f"{prefix}onopy"]:
        del settings[f"{prefix}onopy"]["serial"]

    ctx = PhonopyContext(settings=settings)
    ctx.settings.atoms = atoms.copy()
    if calc is not None and calc.name.lower() != "aims":
        ctx.settings.atoms.set_calculator(calc)
    else:
        aims_ctx = AimsContext(settings=ctx.settings, workdir=ctx.workdir)
        # set reference structure for aims calculation and make sure forces are computed
        aims_ctx.ref_atoms = make_supercell(
            atoms, get_3x3_matrix(ctx.settings.obj.supercell_matrix)
        )
        aims_ctx.settings.obj["compute_forces"] = True
        calc = setup_aims(aims_ctx, False, False)
        ctx.settings.atoms.set_calculator(calc)

    # outputs = bootstrap(name=f"{prefix}onopy", settings=settings, **kwargs_boot)
    outputs = bootstrap(ctx)

    outputs["metadata"]["supercell"] = {"atoms": outputs["metadata"]["atoms"]}
    outputs["metadata"]["primitive"] = input2dict(atoms)
    outputs["prefix"] = prefix
    outputs["settings"] = ph_settings.copy()
    return outputs


def bootstrap_phonon(
    atoms, calc, kpt_density=None, ph_settings=None, ph3_settings=None, fw_settings=None
):
    """Creates a Settings object and passes it to the bootstrap function

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Atoms object of the primitive cell
    calc: ase.calculators.calulator.Calculator
        Calculator for the force calculations
    kpt_density: float
        k-point density for the MP-Grid
    ph_settings: dict
        kwargs for phonopy setup
    ph3_settings: dict
        kwargs for phono3py setup
    fw_settings: dict
        FireWork specific settings

    Returns
    -------
    outputs: dict
        The output of vibes.phonopy.workflow.bootstrap for phonopy and phono3py
    """
    settings = TaskSettings(name=None, settings=Settings(settings_file=None))
    settings.atoms = atoms

    if kpt_density:
        settings["control_kpt"] = AttributeDict({"density": kpt_density})

    outputs = []
    at = atoms.copy()
    at.set_calculator(None)

    if ph_settings:
        outputs.append(setup_phonon_outputs(ph_settings, settings, "ph", at, calc))

    if ph3_settings:
        outputs.append(setup_phonon_outputs(ph3_settings, settings, "ph3", at, calc))

    if kpt_density:
        for out in outputs:
            out["kpt_density"] = kpt_density
    return outputs


def collect_to_trajectory(workdir, trajectory, calculated_atoms, metadata):
    """Collects forces to a single trajectory file

    Parameters
    ----------
        workdir: str
            working directory for the task
        trajectory: str
            file name for the trajectory file
        calculated_atoms: list of dicts
            Results of the force calculations
        metadata: dict
            metadata for the phonon calculations
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
        temp_atoms = []
        for atoms_dict in calculated_atoms:
            calc_dict = atoms_dict.pop("calculator_dict")
            temp_atoms.append(dict2atoms(atoms_dict, calc_dict))
    else:
        temp_atoms = calculated_atoms.copy()
    try:
        calculated_atoms = sorted(
            temp_atoms,
            key=lambda x: x.info[displacement_id_str]
            if x
            else len(calculated_atoms) + 1,
        )
    except KeyError:
        calculated_atoms = sorted(
            temp_atoms,
            key=lambda x: int(x.info["info_str"][1].split("T = ")[1].split(" K")[0])
            if x
            else len(calculated_atoms) + 1,
        )
    for atoms in calculated_atoms:
        if atoms:
            step2file(atoms, atoms.calc, str(traj))


def phonon_postprocess(func_path, phonon_times, kpt_density, **kwargs):
    """Performs phonon postprocessing steps

    Parameters
    ----------
    func_path: str
        The path to the postprocessing function
    phonon_times: list of ints
        The time it took to calculate the phonon forces in seconds
    kwargs: dict
        The keyword arguments for the phonon calculations

    Returns
    -------
    phonopy.Phonopy or phono3py.phonon3.Phono3py
        The Phonopy or Phono3py object generated by the post processing
    """
    func = get_func(func_path)
    return func(**kwargs)


def prepare_gruneisen(settings, primitive, vol_factor):
    """Prepare a Gruneisen calculation

    Parameters
    ----------
    settings: Settings
        The settings object for the calculation
    primitive: ase.atoms.Atoms
        The primitive cell for the phonon calculation
    vol_factor: float
        The volume rescaling factor

    Returns
    -------
    Workflow:
        A Fireworks workflow for the gruneisen calculation
    """
    dist_primitive = primitive.copy()
    scaled_pos = dist_primitive.get_scaled_positions()
    dist_primitive.cell *= vol_factor ** (1.0 / 3.0)
    dist_primitive.set_scaled_positions(scaled_pos)

    for constraint in dist_primitive.constraints:
        if isinstance(constraint, FixCartesianParametricRelations):
            constraint.params = []
            constraint.Jacobian = np.zeros((constraint.Jacobian.shape[0], 0))
            constraint.Jacobian_inv = np.zeros((0, constraint.Jacobian.shape[0]))
            constraint.const_shift = dist_primitive.cell.flatten()

    dist_settings = Settings()
    for sec_key, sec_val in settings.items():
        if isinstance(sec_val, dict):
            dist_settings[sec_key] = AttributeDict()
            for key, val in sec_val.items():
                try:
                    dist_settings[sec_key][key] = val.copy()
                except AttributeError:
                    dist_settings[sec_key][key] = val
        else:
            dist_settings[sec_key] = val

    file_original = dist_settings.geometry.pop("file", None)
    dist_primitive.write(
        "geometry.in.temp", format="aims", geo_constrain=True, scaled=True
    )
    dist_settings.geometry["file"] = "geometry.in.temp"
    dist_primitive.calc = setup_aims(
        ctx=AimsContext(settings=dist_settings), verbose=False, make_species_dir=False
    )
    Path("geometry.in.temp").unlink()
    dist_settings.geometry["file"] = file_original

    dist_settings.atoms = dist_primitive

    return generate_workflow(dist_settings, dist_primitive, launchpad_yaml=None)


def setup_gruneisen(settings, trajectory, constraints, _queueadapter, kpt_density):
    """Set up the finite difference gruniesen parameter calculations

    Parameters
    ----------
    settings: dict
        The workflow settings
    trajectory: str
        The trajectory file for the equilibrium phonon calculation
    constraints: list of dict
        list of relevant constraint dictionaries for relaxations
    _queueadapter: dict
        The queueadapter for the equilibrium phonon calculations
    kpt_density:
        The k-point density to use for the calculations

    Returns
    -------
    pl_gruneisen: Workflow
        Workflow to run the positive volume change Gruniesen parameter calcs
    mn_gruneisen: Workflow
        Workflow to run the positive volume change Gruniesen parameter calcs
    """
    # Prepare settings by reset general work_dir and do not reoptimize k_grid
    settings["general"]["opt_kgrid"] = False
    settings["phonopy"]["get_gruniesen"] = False
    settings["phonopy"]["converge_phonons"] = False

    settings.pop("statistical_sampling", None)

    if _queueadapter:
        settings["phonopy_qadapter"] = _queueadapter

    # Get equilibrium phonon
    eq_phonon = postprocess(trajectory)
    _, metadata = reader(trajectory, get_metadata=True)

    settings["phonopy"]["supercell_matrix"] = eq_phonon.get_supercell_matrix()
    settings["phonopy"]["symprec"] = metadata["Phonopy"].get("symprec", 1e-5)
    settings["phonopy"]["displacement"] = metadata["Phonopy"]["displacements"][0][1]

    settings["control"] = {}
    for key, val in metadata["calculator"]["calculator_parameters"].items():
        settings["control"][key] = val
    settings["control"].pop("kgrid", None)
    settings["control_kpt"] = {"density": kpt_density}

    primitive = to_Atoms(eq_phonon.get_primitive())
    add_constraints = []
    for constr in constraints:
        constraint = dict2constraint(constr)
        if isinstance(constraint, FixScaledParametricRelations):
            if not constraint.params:
                settings["general"]["relax_structure"] = False
                settings["control"].pop("relax_unit_cell", None)
                settings.pop("relaxation", None)
            else:
                add_constraints.append(constraint)

    primitive.constraints = add_constraints

    settings["general"]["relax_structure"] = True
    if "relaxation" not in settings:
        settings["relaxation"] = {}
    settings["control"]["relax_unit_cell"] = "none"
    settings["relaxation"]["relax_unit_cell"] = "none"

    pl_gruneisen = prepare_gruneisen(settings, primitive, 1.01)
    mn_gruneisen = prepare_gruneisen(settings, primitive, 0.99)

    return pl_gruneisen, mn_gruneisen
