""" Provide an aims calculator without much ado """
from os import path
from pathlib import Path

# from ase.calculators.aims import Aims
from hilde.calculators.aims_calc import Aims
from hilde.settings import Settings


def setup_aims(
    custom_settings={}, workdir=None, settings=None, config_file="../../hilde.cfg"
):
    """Set up an aims calculator.

    Args:
        custom_settings (dict): Settings that replace the minimal defaults
        workdir (str): directory to work in
        config_file (str): path to config file

    Returns:
        Aims: ASE calculator object

    """

    if settings is None:
        settings = Settings(config_file)

    default_settings = {**settings.control}
    ase_settings = {
        "aims_command": settings.machine.aims_command,
        "species_dir": path.join(settings.machine.basissetloc, settings.basisset.type),
    }

    # Check if basisset type is supposed to be changed by custom settings
    if "species_type" in custom_settings:
        species_type = custom_settings["species_type"]
        if "species_dir" not in custom_settings:
            species_dir = path.join(settings.machine.basissetloc, species_type)
        else:
            species_dir = path.join(custom_settings["species_dir"], species_type)

        custom_settings["species_dir"] = species_dir
        del (custom_settings["species_type"])

    if "socketio" in settings and settings.socketio.port is not None:
        custom_settings.update(
            {"use_pimd_wrapper": ("localhost", settings.socketio.port)}
        )
        if "use_socketio" in custom_settings:
            del custom_settings["use_socketio"]

    aims_settings = {**default_settings, **custom_settings, **ase_settings}

    if workdir:
        return Aims(label=Path(workdir).absolute(), **aims_settings)
    else:
        return Aims(**aims_settings)
