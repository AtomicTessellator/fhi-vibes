"""Wrappers to prepare statistical sampling"""

from hilde.fireworks.tasks.phonopy_phono3py_functions import setup_calc
from hilde.helpers.attribute_dict import AttributeDict
from hilde.phonon_db.ase_converters import calc2dict
from hilde.settings import Settings
from hilde.statistical_sampling.workflow import bootstrap as bootstrap_stat_samp


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
        (dict): The output of hilde.statistical_sampling.workflow.bootstrap
    """
    settings = Settings(settings_file=None)
    settings.atoms = atoms
    if kpt_density:
        settings["control_kpt"] = AttributeDict({"density": kpt_density})

    outputs = []

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
            atoms=atoms, name="statistical_sampling", settings=settings, **kwargs_boot
        )
        stat_samp_out["prefix"] = "stat_samp"
        stat_samp_out["settings"] = stat_samp_settings.copy()
        stat_samp_out["calculator"] = calc2dict(calc)
        outputs.append(stat_samp_out)
    else:
        raise IOError("Stat sampling requires a settings object")
    return outputs
