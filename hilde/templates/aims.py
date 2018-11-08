""" Provide an aims calculator without much ado """

from pathlib import Path
from ase.calculators.aims import Aims
from hilde.settings import Settings


def setup_aims(
    port=None, custom_settings={}, workdir=None, config_file="../../hilde.cfg"
):
    """Set up an aims calculator.

    Args:
        custom_settings (dict): Settings that replace the minimal defaults
        workdir (str): directory to work in
        config_file (str): path to config file

    Returns:
        Aims: ASE calculator object

    """

    settings = Settings(config_file)
    # Check if basisset type is supposed to be changed by custom settings
    if "species_type" in custom_settings:
        species_type = custom_settings["species_type"]
        species_dir = str(Path(settings.machine.basissetloc) / species_type)
        custom_settings["species_dir"] = species_dir
        del (custom_settings["species_type"])

    default_settings = settings.ase_settings

    if port is not None:
        custom_settings.update({"use_pimd_wrapper": ("localhost", port)})

    aims_settings = {**default_settings, **custom_settings}

    if workdir:
        return Aims(label=Path(workdir).absolute(), **aims_settings)
    else:
        return Aims(**aims_settings)
