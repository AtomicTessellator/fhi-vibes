"""Generators for Jobs for the balsam generators"""
from ase.calculators.aims import Aims

from balsam.launcher.dag import BalsamJob

import numpy as np

from hilde.aims.context import AimsContext

from hilde.balsam.data_encoder import encode

from hilde.helpers.hash import hash_atoms
from hilde.phonon_db.ase_converters import atoms2dict, calc2dict
from hilde.phonopy.context import PhonopyContext
from hilde.phonopy import defaults as ph_defaults


def generate_aims_job(settings, basisset=None, calc_number=0):
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
    ctx = AimsContext(settings=settings, read_config=False)

    for sec_key in ["phonopy", "phono3py", "md"]:
        if sec_key in ctx.settings:
            del (ctx.settings[sec_key])

    atoms = settings.atoms.copy()
    calc = Aims(**settings.atoms.calc.parameters)
    atoms.set_calculator(None)

    if "relaxation" in settings:
        method = settings.relaxation.get("method", "trm")
        force_crit = settings.relaxation.get("conv_crit", 1e-3)
        relax_unit_cell = settings.relaxation.get("relax_unit_cell", "full")
    else:
        method = "trm"
        force_crit = 1e-3
        relax_unit_cell = "full"

    if basisset:
        ctx.settings["basisset"]["type"] = basisset
    else:
        basisset = ctx.settings["basisset"]["type"]

    update_control = {
        "relax_geometry": f"{method} {force_crit}",
        "scaled": True,
        "relax_unit_cell": relax_unit_cell,
        "use_sym": True,
    }

    calc.parameters.update(update_control)
    ctx.settings.control.update(update_control)

    data = {
        "atoms": atoms2dict(atoms),
        "calculator": calc2dict(calc),
        "ctx": ctx.settings,
        "calc_number": calc_number,
    }

    if f"{basisset}_rel_qadapter" in settings:
        num_nodes = settings[f"{basisset}_rel_qadapter"]["nodes"]
    else:
        num_nodes = 1

    job_suffix = f"{atoms.symbols.get_chemical_formula()}/{hash_atoms(atoms)}"
    return BalsamJob(
        name=f"aims_calc_{basisset}_0",
        workflow=f"{job_suffix}",
        description=f"Aims calculation for {atoms.symbols.get_chemical_formula()} with hash {hash_atoms(atoms)}.",
        application="hilde-run-aims",
        num_nodes=num_nodes,
        ranks_per_node=1,
        cpu_affinity="depth",
        data=encode(data),
        input_files="",
    )


def generate_phonopy_job(settings):
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
        if sec_key in ctx.settings:
            del (ctx.settings[sec_key])

    atoms = settings.atoms.copy()
    calc = Aims(**settings.atoms.calc.parameters)
    atoms.set_calculator(None)

    data = {
        "atoms": atoms2dict(atoms),
        "calculator": calc2dict(calc),
        "ctx": ctx.settings,
        "prev_dos_fp": None,
        "sc_matrix_original": list(ctx.settings.phonopy.supercell_matrix.flatten()),
        "conv_crit": ctx.settings.phonopy.conv_crit,
    }
    basisset = ctx.settings.basisset.type
    if f"{basisset}_rel_qadapter" in settings:
        num_nodes = settings[f"{basisset}_rel_qadapter"]["nodes"]
    else:
        num_nodes = 1

    natoms = int(
        round(len(atoms) * np.linalg.det(ctx.settings.phonopy.supercell_matrix))
    )

    job_suffix = f"{atoms.symbols.get_chemical_formula()}/{hash_atoms(atoms)}"
    return BalsamJob(
        name=f"phonopy_{basisset}_natoms_{natoms}",
        workflow=f"{job_suffix}",
        description=f"phonopy calculation for {atoms.symbols.get_chemical_formula()} with hash {hash_atoms(atoms)} and {natoms} supercell.",
        application="hilde-run-phonopy",
        num_nodes=num_nodes,
        ranks_per_node=1,
        cpu_affinity="depth",
        data=encode(data),
        input_files="",
    )


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
        del (sym_params[1:n_lat + 1])
        symmetry_block[1] = " ".join(sym_params)

    jobs = []
    if "basisset" in settings.general:
        basis = settings.general.pop("basisset")
    else:
        basis = "light"

    settings.atoms = atoms

    if not symmetry_block or int(sym_nparams[1]) - n_lat > 0:
        jobs.append(generate_aims_job(settings, "light"))
        # Tighter Basis Set Relaxation
        use_tight_relax = settings.general.get("use_tight_relax", False)
        if basis != "light" or use_tight_relax:
            if use_tight_relax:
                basisset_type = "tight"
            else:
                basisset_type = basis
            jobs.append(generate_aims_job(settings, basisset_type))

    ignore_keys = ["trigonal", "q_mesh"]
    for key, val in ph_defaults.items():
        if key not in settings.phonopy and key not in ignore_keys:
            settings.phonopy[key] = val

    settings.phonopy["basisset_type"] = basis

    if "supercell_matrix" not in settings.phonopy:
        raise IOError("Initial supercell_matrix must be provided")

    settings.phonopy.converge_phonons = False
    settings.phonopy.get_gruneisen = False

    jobs.append(generate_phonopy_job(settings))

    job_dep = {}
    for ii in range(len(jobs) - 1):
        job_dep[jobs[ii]] = jobs[ii + 1]
    return jobs, job_dep
