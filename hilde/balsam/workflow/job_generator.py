"""Generators for Jobs for the balsam generators"""
from pathlib import Path

import ase
from ase.calculators.aims import Aims
from ase import units as u
from ase.md.velocitydistribution import PhononHarmonics

from balsam.core.models import BalsamJob

import balsam.launcher.dag as dag

import numpy as np

from hilde.aims.context import AimsContext
from hilde.aims.setup import setup_aims

from hilde.balsam.data_encoder import decode, encode
from hilde import konstanten as const
from hilde.helpers.converters import input2dict
from hilde.helpers.hash import hash_atoms
from hilde.helpers.numerics import get_3x3_matrix
from hilde.helpers.supercell import make_supercell
from hilde.helpers.watchdogs import str2time
from hilde.phonon_db.ase_converters import atoms2dict, calc2dict
from hilde.phonopy.context import PhonopyContext
from hilde.phonopy.utils import remap_force_constants
from hilde.phonopy.workflow import bootstrap
from hilde.settings import Settings
from hilde.structure.convert import to_Atoms

ranks_per_node = Settings()["comp_node_stats"]["cores_per_node"]
max_nodes = Settings()["comp_node_stats"]["max_nodes"]


def generate_aims_job(settings, basisset=None, calc_number=0, app="run-aims"):
    """Generate an Aims calculation

    Parameters
    ----------
    settings: Settings
        The settings object for the workflow

    Returns
    -------
    BalsamJob:
        Job for aims calculation
    """
    # Setup the context
    ctx = AimsContext(settings=settings, read_config=False)

    for sec_key in ["phonopy", "phono3py", "md"]:
        ctx.settings.pop(sec_key, None)

    atoms = settings.atoms.copy()
    calc = Aims(**settings.atoms.calc.parameters)
    atoms.set_calculator(None)

    # Get the correct basiss settigns
    if basisset:
        ctx.settings["basisset"]["type"] = basisset
    else:
        basisset = ctx.settings["basisset"]["type"]

    # Check if a relaxation is being done, and if so get the relevant parameters
    if "relaxation" in settings:
        settings.general["relax_structure"] = True
        method = settings.relaxation.get("method", "trm")
        force_crit = settings.relaxation.get("conv_crit", 1e-3)
        relax_unit_cell = settings.relaxation.get("relax_unit_cell", "full")
    elif settings.general.get("relax_structure", True):
        method = "trm"
        force_crit = 1e-3
        relax_unit_cell = "full"

    # Update the control settings in the context and calculator
    update_control = dict()
    if settings.general.get("relax_structure", True):
        update_control["relax_geometry"] = f"{method} {force_crit}"
        update_control["relax_unit_cell"] = relax_unit_cell
    calc.parameters.update(update_control)
    ctx.settings.control.update(update_control)

    # Generate the initial data for teh job
    data = {
        "atoms": atoms2dict(atoms),
        "calculator": calc2dict(calc),
        "ctx": ctx.settings,
        "calc_number": calc_number,
    }
    if settings.general.get("relax_structure", True):
        data["relaxation"] = {
            "method": method,
            "force_crit": force_crit,
            "relax_unit_cell": relax_unit_cell,
        }

    # Get initial job parameters
    natoms = len(atoms)
    if f"{basisset}_rel_qadapter" in settings:
        num_nodes = settings[f"{basisset}_rel_qadapter"].get("nodes", 1)
        walltime = settings[f"{basisset}_rel_qadapter"].get("walltime", "1:00:00")
        walltime = str2time(walltime) / 60
    else:
        num_nodes = int(np.ceil(1.0 * natoms / ranks_per_node))
        walltime = 30

    # Create and save the job
    job_suffix = f"{atoms.symbols.get_chemical_formula()}/{hash_atoms(atoms)}"
    job = BalsamJob(
        name=f"aims_calc_{basisset}_0",
        workflow=f"{job_suffix}",
        description=f"Aims calculation for {atoms.symbols.get_chemical_formula()} with hash {hash_atoms(atoms)}.",
        application=app,
        num_nodes=min(num_nodes, max_nodes),
        wall_time_minutes=walltime,
        ranks_per_node=ranks_per_node,
        cpu_affinity="depth",
        data=encode(data),
        input_files="",
    )
    job.save()
    return job


def generate_phonopy_jobs(settings, data):
    """Generate a list of phonopy supercell calculation

    Parameters
    ----------
    settings: PhonopySettings
        The context of the phonopy calculation
    data: dict
        The decoded dag.current_job.data

    Returns
    -------
    list(BalsamJobs):
        The list of phonopy calculations
    """
    # Get the initial context
    ctx = PhonopyContext(settings=settings, read_config=False)

    # get the calculator and atoms object from the ctx
    calc = setup_aims(AimsContext(settings, read_config=True))
    atoms = settings.atoms.copy()

    ctx.ref_atoms = atoms
    ph_initialization = bootstrap(ctx)

    chem_formula = atoms.symbols.get_chemical_formula()
    atoms_hash = hash_atoms(atoms)

    workflow = f"{chem_formula}/{atoms_hash}"

    atoms.set_calculator(calc)

    natoms = int(
        round(len(atoms) * np.linalg.det(ctx.settings.phonopy.supercell_matrix))
    )

    copy_wd_base = f"{settings.general.workdir_cluster}/{workflow}/sc_natoms_{natoms}/"
    name = f"aims_calc_ph_{natoms}_supercell_"
    description = f"FHI-aims calculation for the forces of a {natoms} atom supercell of {chem_formula} and hash {atoms_hash}. Displacement number: "

    ph_update_data = {
        "metadata": ph_initialization["metadata"],
        "number_of_sc": len(ph_initialization["atoms_to_calculate"]),
        "analysis_wd": copy_wd_base,
    }
    data["ph_data"].update(ph_update_data)

    Path(copy_wd_base).mkdir(exist_ok=True, parents=True)
    ctx.settings.write(f"{copy_wd_base}/phonopy.in")
    atoms.write(
        f"{copy_wd_base}/geometry.in", format="aims", scaled=True, geo_constrain=True
    )
    job_list = []

    for ii, atoms in enumerate(ph_initialization["atoms_to_calculate"]):
        if atoms is None:
            continue
        ctx = AimsContext(settings=settings, read_config=False)
        ctx.settings.atoms = atoms
        ctx.settings.atoms.set_calculator(calc)

        job = generate_aims_job(ctx.settings, app="run-aims-post-phonopy")
        job.num_nodes = min(data["ph_data"]["nodes"], max_nodes)
        job.wall_time_minutes = data["ph_data"]["walltime"]
        job.workflow = workflow
        job.name = name + f"{ii}"
        job.description = description + f"{ii}"

        job_data = decode(job.data)
        job_data["ph_data"] = data["ph_data"].copy()
        job_data["analysis_wd"] = ph_update_data["analysis_wd"] + f"{ii:05d}"
        if "relaxation" in data:
            job_data["relaxation"] = data["relaxation"]
        job.data = encode(job_data)

        job.save()

        job_list.append(job)
    return job_list


def get_phonon_setup_data(settings, qadapter=None, prev_dos_fp=None):
    """Get the ph_data needed to setup a phonopy calculation

    Parameters
    ----------
    settings: Settings
        Settings for the phonopy task
    qadapter: dict
        Job runtime properties
    prev_dos_fp: hilde.materials_fp.material_fingerprint.fp_tup
        Current Phonon DOS fingerprint (previous for next step)

    Returns
    -------
    ph_data: dict
        The needed information to postprocess the phonopy step
    """

    atoms = settings.atoms
    sc_matrix = get_3x3_matrix(settings.phonopy.get("supercell_matrix"))
    natoms = len(atoms) * np.linalg.det(sc_matrix)

    num_nodes = int(np.ceil(1.0 * natoms / ranks_per_node))
    walltime = "00:30:00"

    if qadapter is not None:
        num_nodes = qadapter.get("nodes", num_nodes)
        walltime = qadapter.get("walltime", walltime)
    walltime = int(np.ceil(float(str2time(walltime)) / float(60)))

    ph_data = {
        "ph_primitive": atoms2dict(atoms),
        "prev_dos_fp": prev_dos_fp,
        "walltime": walltime,
        "nodes": num_nodes,
    }
    return ph_data


def make_job_phonopy_setup(settings, job):
    """Generate an Aims calculation

    Parameters
    ----------
    settings: Settings
        The settings object for the workflow

    Returns
    -------
    BalsamJob:
        Job for aims calculation
    """
    ctx = PhonopyContext(settings=settings, read_config=False)

    for sec_key in ["phono3py", "md"]:
        ctx.settings.pop(sec_key, None)

    ctx.settings.atoms = settings.atoms.copy()
    ctx.settings.atoms.set_calculator(None)

    ph_data = get_phonon_setup_data(
        ctx.settings, ctx.settings.pop("phonopy_qadapter", None)
    )
    sc_matrix = settings.phonopy.get("supercell_matrix")
    del ph_data["ph_primitive"]
    ph_data["sc_matrix_original"] = sc_matrix
    ph_data["conv_crit"] = settings.phonopy.get("conv_crit", 0.95)
    ph_data["ctx"] = ctx.settings

    job.application = "run-aims-setup-phonopy"
    job_data = decode(job.data)
    job_data["ph_data"] = ph_data
    job.data = encode(job_data)
    job.save()

    return job


def generate_gruneisen_jobs(settings, vol_factor):
    """Generate a FireWork to calculate the Gruniesen Parameter with finite differences

    Parameters
    ----------
    sdettings: Settings
        The Workflow Settings

    Returns
    -------
    BalsamJob:
        The Gruniesen
    """

    atoms = settings.atoms.copy()
    atoms.set_calculator(settings.atoms.calc)

    atoms.cell *= vol_factor
    atoms.positions *= vol_factor

    symmetry_block = atoms.info.get("symmetry_block", None)

    if symmetry_block:
        sym_nparams = symmetry_block[0].split()
        n_lat = int(sym_nparams[2])
        for ii in range(1, 3):
            sym_nparams[ii] = str(int(sym_nparams[ii]) - n_lat)
        symmetry_block[0] = " ".join(sym_nparams)
        sym_params = symmetry_block[1].split()
        del sym_params[1 : n_lat + 1]
        symmetry_block[1] = " ".join(sym_params)
        for ii in range(3):
            symmetry_block[
                ii + 2
            ] = f"symmetry_lv {atoms.cell[ii, 0]}, {atoms.cell[ii, 1]}, {atoms.cell[ii, 2]}"
        atoms.info["symmetry_block"] = symmetry_block
    jobs = []
    if "basisset" in settings.general:
        basis = settings.general.pop("basisset")
    else:
        basis = "light"

    settings.atoms = atoms
    relaxation = decode(dag.current_job.data).get("relaxation", None)
    if (not symmetry_block or int(sym_nparams[-1]) > 0) and relaxation is not None:
        settings.general["relax_structure"] = True
        settings["relaxation"] = relaxation
        jobs.append(generate_aims_job(settings, "light"))
        # Tighter Basis Set Relaxation
        use_tight_relax = settings.general.get("use_tight_relax", False)
        if basis != "light" or use_tight_relax:
            if use_tight_relax:
                basisset_type = "tight"
            else:
                basisset_type = basis
            jobs.append(generate_aims_job(settings, basisset_type))
        make_job_phonopy_setup(settings, jobs[-1])
        return jobs

    ph_data = get_phonon_setup_data(settings, settings.pop("phonopy_qadapter", None))
    ph_data["ctx"] = settings
    data = dict(ph_data=ph_data)
    return generate_phonopy_jobs(settings, data)


def generate_stat_samp_jobs(settings, phonon):
    """Generate a FireWork to calculate the Gruniesen Parameter with finite differences

    Parameters
    ----------
    sdettings: Settings
        The Workflow Settings

    Returns
    -------
    BalsamJob:
        The Gruniesen
    """
    atoms = to_Atoms(phonon.get_primitive())
    atoms.set_calculator(settings.atoms.calc)

    if "supercell_matrix" in settings.statistical_sampling:
        sc = make_supercell(atoms, settings.statistical_sampling.supercell_matrix)
    else:
        sc = to_Atoms(phonon.get_supercell())

    force_constants = phonon.get_force_constants().copy()

    temperatures = settings.statistical_sampling.get("temperatures", list())
    # If using Debye temperature calculate it
    debye_temp_fact = settings.statistical_sampling.get("debye_temp_fact", list())
    if debye_temp_fact:
        if phonon is None:
            raise IOError(
                "Debye Temperature must be calculated with phonopy, please add phonon_file"
            )
        phonon.set_mesh([51, 51, 51])
        phonon.set_total_DOS(freq_pitch=0.01)
        phonon.set_Debye_frequency()
        debye_temp = phonon.get_Debye_frequency() * const.THzToEv / const.kB
        temperatures += [tt * debye_temp for tt in debye_temp_fact]

    if not temperatures:
        raise IOError("temperatures must be given to do harmonic analysis")

    if "deterministic" in settings.statistical_sampling:
        deterministic = settings.statistical_sampling["deterministic"]

    nsamples = settings.statistical_sampling.get("n_samples", 1)
    # Generate metadata
    metadata = {
        "ase_vesrsion": ase.__version__,
        "n_samples_per_temperature": nsamples,
        "temperatures": temperatures,
        "deterministic": deterministic,
        "force_constants": force_constants,
        "supercell": input2dict(sc),
        "primitive": input2dict(atoms),
        **input2dict(sc),
    }

    rng_seed = settings.statistical_sampling.get(
        "rng_seed", np.random.randint(2 ** 32 - 1)
    )
    if not deterministic:
        if not isinstance(rng_seed, int):
            rng_seed = int(rng_seed)
        metadata["rng_seed"] = rng_seed
    rng = np.random.RandomState(rng_seed)

    force_constants = remap_force_constants(
        force_constants, atoms, to_Atoms(phonon.get_supercell()), sc, two_dim=True
    )
    assert force_constants.shape[0] == force_constants.shape[1] == len(sc) * 3

    job_list = []
    # Set up thermally displaced supercells
    td_cells = list()
    for temp in temperatures:
        sc.info["temperature"] = temp
        td_cells += prepare_phonon_harmonic_sampling(
            sc, force_constants, temp, nsamples, deterministic, rng
        )
    for cell in td_cells:
        ctx = AimsContext(settings=settings, read_config=False)
        ctx.settings.atoms = cell.copy()
        ctx.settings.atoms.calc = setup_aims(AimsContext(settings=settings))
        job = generate_aims_job(ctx.settings)
        job.workflow = dag.current_job.workflow
        job.save()
        job_list.append(job)

    return job_list


def prepare_phonon_harmonic_sampling(
    atoms, force_constants, temperature, n_samples=1, deterministic=True, rng=np.random
):
    """
    Generates a list of displaced supercells based on a thermal excitation of phonons
    Parameters:
        atoms (ase.atoms.Atoms): Non-displaced supercell
        force_constants(np.ndarray(shape=(3*len(atoms), 3*len(atoms)))): Force constant matrix for atoms
        temperature (float): Temperature to populate the phonons modes at
        n_samples(int): number of samples to generate
        deterministic(bool): If True then displace atoms with +/- the amplitude according to PRB 94, 075125
    Returns (list of ase.atoms.Atoms): The thermally displaced supercells
    """
    thermally_disp_cells = []
    for _ in range(n_samples):
        td_cell = atoms.copy()
        PhononHarmonics(
            td_cell,
            force_constants,
            quantum=False,
            temp=temperature * u.kB,
            rng=rng,
            plus_minus=deterministic,
            failfast=True,
        )
        thermally_disp_cells.append(td_cell)
    return thermally_disp_cells
