"""Wrappers to prepare statistical sampling"""
import ase

import numpy as np

from vibes.helpers.attribute_dict import AttributeDict
from vibes.helpers.converters import input2dict, dict2atoms, calc2dict
from vibes.helpers.k_grid import k2d, update_k_grid
from vibes.helpers.supercell import make_supercell
from vibes.helpers.warnings import warn
from vibes.phonopy.postprocess import postprocess as postprocess_ph
from vibes.phonopy.utils import get_force_constants_from_trajectory, remap_force_constants
from vibes.phonopy.wrapper import preprocess as get_debye_temperature
from vibes.scripts.create_samples import generate_samples
from vibes.settings import TaskSettings, Settings
from vibes.structure.convert import to_Atoms
from vibes.trajectory import reader


def bootstrap(atoms, name="statistical_sampling", settings=None, **kwargs):
    """ load settings, prepare atoms, calculator, and phonopy """

    if settings is None:
        settings = TaskSettings(name=None, settings=Settings())

    stat_sample_settings = {}

    if name not in settings:
        warn(f"Settings do not contain {name} instructions.", level=1)
    else:
        stat_sample_settings.update(settings[name])

    _, ph_metadata = reader(stat_sample_settings["phonon_file"], get_metadata=True)

    atoms_dict = ph_metadata["atoms"].copy()
    ph_atoms = dict2atoms(atoms_dict, ph_metadata["calculator"], False)
    calc = ph_atoms.calc
    kpt_density = k2d(ph_atoms, calc.parameters["k_grid"])

    # Get sampling metadata
    stat_sample_settings.update(kwargs)
    metadata = get_metadata(**stat_sample_settings)

    # Generate Samples
    sc = dict2atoms(metadata["supercell"])
    fc = remap_force_constants(
        metadata["force_constants"],
        dict2atoms(metadata["primitive"]),
        sc,
        sc,
        two_dim=True,
    )

    td_cells = []
    for temp in metadata["temperatures"]:
        td_cells += generate_samples(
            sc, temp, force_constants=fc, **metadata["generate_sample_args"]
        )

    calc = update_k_grid(td_cells[0], calc, kpt_density)
    calc.parameters.pop("scaled", False)
    # save metadata
    metadata["calculator"] = calc2dict(calc)

    return {
        "atoms_to_calculate": td_cells,
        "calculator": calc,
        "metadata": metadata,
        "workdir": name,
        "settings": stat_sample_settings,
    }


def get_metadata(phonon_file, temperatures=None, debye_temp_fact=None, **kwargs):
    """
    Sets ups the statistical sampling
    Parameters:
        phonon_file (str): String to the phonopy trajectory
        temperatures (list (floats)): List of temperatures to excite the phonons to
        debye_temp_fact (list(floats)): List of factors to multiply the debye temperature by to populate temperatures
    Additional arguments in kwargs:
        supercell_matrix (np.ndarray(int)): Supercell matrix used to create the supercell from atoms
        n_samples (int): number of samples to create (default: 1)
        mc_rattle (bool): hiphive mc rattle
        quantum (bool): use Bose-Einstein distribution instead of Maxwell-Boltzmann
        deterministic (bool): create sample deterministically
        sobol (bool): Use sobol numbers for the sampling
        rng_seed (int): seed for random number generator
    Returns: (list(dicts)): A list of thermally displaced supercells to calculate the forces on
    """
    if temperatures is None:
        temperatures = ()

    # Set up supercell and Force Constants
    phonon = postprocess_ph(
        phonon_file, write_files=False, calculate_full_force_constants=False
    )
    atoms = to_Atoms(phonon.get_primitive())
    if "supercell_matrix" in kwargs:
        sc = make_supercell(atoms, np.array(kwargs["supercell_matrix"]).reshape(3, 3))
    else:
        sc = to_Atoms(phonon.get_supercell())

    force_constants = get_force_constants_from_trajectory(phonon_file, sc, reduce_fc=True)
    # If using Debye temperature calculate it
    if debye_temp_fact is not None:
        if phonon is None:
            raise IOError(
                "Debye Temperature must be calculated with phonopy, please add phonon_file"
            )

        debye_temp = get_debye_temperature(phonon, 5e-3)[0]
        temperatures += [tt * debye_temp for tt in debye_temp_fact]
    elif temperatures is None:
        raise IOError("temperatures must be given to do harmonic analysis")

    # Generate metadata
    metadata = {
        "ase_vesrsion": ase.__version__,
        "n_samples_per_temperature": kwargs.get("n_samples", 1),
        "temperatures": temperatures,
        "force_constants": force_constants,
        "supercell": input2dict(sc)["atoms"],
        "primitive": input2dict(atoms)["atoms"],
        "atoms": input2dict(sc)["atoms"],
    }

    if not kwargs.get("deterministic", True):
        if kwargs.get("rng_seed", None) is None:
            rng_seed = np.random.randint(2 ** 32 - 1)
        else:
            rng_seed = int(kwargs.get("rng_seed"))
    else:
        rng_seed = None

    metadata["generate_sample_args"] = {
        "n_samples": kwargs.get("n_samples", 1),
        "rattle": kwargs.get("mc_rattle", False),
        "quantum": kwargs.get("quantum", False),
        "deterministic": kwargs.get("deterministic", True),
        "plus_minus": kwargs.get("plus_minus", True),
        "gauge_eigenvectors": kwargs.get("gauge_eigenvectors", True),
        "ignore_negative": kwargs.get("ignore_negative", True),
        "sobol": kwargs.get("sobol", False),
        "random_seed": rng_seed,
        "propagate": False,
        "format": "aims",
    }

    return metadata


def bootstrap_stat_sample(
    atoms, calc, kpt_density=None, stat_samp_settings=None, fw_settings=None
):
    """
    Initializes the statistical sampling task
    Args:
        atoms (ASE Atoms Object): Atoms object of the primitive cell
        calc (ASE Calculator): Calculator for the force calculations
        kpt_density (float): k-point density for the MP-Grid
        stat_samp_settings (dict): kwargs for statistical sampling setup
        fw_settings (dict): FireWork specific settings

    Returns:
        (dict): The output of vibes.statistical_sampling.workflow.bootstrap
    """
    settings = Settings(settings_file=None)
    settings.atoms = atoms
    if kpt_density:
        settings["control_kpt"] = AttributeDict({"density": kpt_density})

    outputs = []

    settings["statistical_sampling"] = stat_samp_settings.copy()
    stat_samp_out = bootstrap(atoms=atoms, name="statistical_sampling", settings=settings)
    stat_samp_out["prefix"] = "stat_samp"
    stat_samp_out["settings"] = stat_samp_settings.copy()

    outputs.append(stat_samp_out)
    return outputs


def postprocess_statistical_sampling(**kwargs):
    """Dummy function for post processing of statistical sampling"""
    return None
