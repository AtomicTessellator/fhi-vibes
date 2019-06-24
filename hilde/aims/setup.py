""" Provide an aims calculator without much ado """
import shutil
from pathlib import Path

# from ase.calculators.aims import Aims
from ase.calculators.aims import Aims
from hilde import DEFAULT_CONFIG_FILE
from hilde.helpers import talk
from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.warnings import warn


def create_species_dir(ctx, folder="basissets"):
    """ create a custom bassiset folder for the computation

    Args:
        ctx (AimsContext): aims context
        older (str/Path): folder to store the basisset

    """

    loc = ctx.basisset_location
    settings = ctx.settings

    # if old section with `basisset.type` is used:
    if "basisset" in settings and "type" in settings.basisset:
        default = settings.basisset.type
        return str(loc / default)
    elif "basissets" in settings and "default" in settings.basissets:
        default = settings.basissets.default
    else:
        warn("basissets not specified in settings.file.", level=2)

    # return default if no atom is given for reference
    ref_atoms = ctx.ref_atoms
    if ref_atoms is None:
        default_path = loc / default
        talk(f"no Atoms object given, return default path {default_path} for basissets")
        return str(default_path)

    folder = ctx.workdir / Path(folder)
    folder.mkdir(exist_ok=True, parents=True)

    symbols = ref_atoms.get_chemical_symbols()
    numbers = ref_atoms.symbols.numbers

    dct = {sym: num for (sym, num) in zip(symbols, numbers)}

    key_vals = (
        (key.capitalize(), val)
        for (key, val) in settings.basissets.items()
        if "default" not in key
    )

    if len(settings.basissets) > 1:
        for (key, val) in key_vals:
            # copy the respective basisset
            shutil.copy(loc / val / f"{dct[key]:02d}_{key}_default", folder)
            del dct[key]

    # add remaining ones
    for key in dct.keys():
        # copy the respective basisset
        shutil.copy(loc / default / f"{dct[key]:02d}_{key}_default", folder)

    return str(folder.absolute())


def setup_aims(ctx):
    """Set up an aims calculator.

    Args:
        ctx (AimsContext): aims context

    Returns:
        Aims: ASE calculator object
    """

    settings = ctx.settings

    aims_settings = settings.obj

    ase_settings = {"aims_command": settings.machine.aims_command}

    if "socketio" in settings and settings.socketio.port is not None:
        aims_settings.update(
            {"use_pimd_wrapper": ("localhost", settings.socketio.port)}
        )

    # create basissetfolder
    species_dir = create_species_dir(ctx)
    ase_settings["species_dir"] = species_dir

    aims_settings = {**aims_settings, **ase_settings}

    calc = Aims(**aims_settings)

    # update k_grid
    if ctx.ref_atoms and "control_kpt" in settings:
        update_k_grid(ctx.ref_atoms, calc, settings.control_kpt.density)

    if "k_grid" not in calc.parameters:
        talk("No k_grid in aims calculator. Check!")

    return calc
