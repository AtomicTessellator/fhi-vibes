""" Provide an aims calculator without much ado """
import shutil
from pathlib import Path

# from ase.calculators.aims import Aims
from ase.calculators.aims import Aims
from hilde.helpers import talk
from hilde.helpers.k_grid import d2k
from hilde.helpers.warnings import warn


choices = ("light", "intermediate", "tight")


class BasissetError(RuntimeError):
    """Raise when the basisset was set up incorrectly"""


def create_species_dir(ctx, folder="basissets", fallback="light"):
    """ create a custom bassiset folder for the computation

    Parameters
    ----------
    ctx: AimsContext
        The context for the calculation
    folder: str or Path
        Folder to store the basisset

    Returns
    -------
    str
        The absolute file path to the species directory
    """

    loc = ctx.basisset_location
    settings = ctx.settings

    # if old section with `basisset.type` is used:
    if "basisset" in settings and "type" in settings.basisset:
        default = settings.basisset.type
        return str(loc / default)
    if "basissets" in settings and "default" in settings.basissets:
        default = settings.basissets.default
        if "fallback" in settings.basissets and settings.basissets.fallback in choices:
            fallback = settings.basissets.fallback
    else:
        warn("basissets not specified in settings.file.", level=2)

    if default not in choices:
        raise BasissetError(f"Species default '{default}' unknown.")

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
        if key not in ("default", "fallback")
    )

    if len(settings.basissets) > 1:
        for (key, val) in key_vals:
            # copy the respective basisset
            add_basisset(loc, val, key, dct[key], folder, fallback=fallback)
            del dct[key]

    # add remaining ones
    for key in dct.keys():
        # copy the respective basisset
        add_basisset(loc, default, key, dct[key], folder, fallback=fallback)

    return str(folder.absolute())


def add_basisset(loc, typ, elem, num, folder, fallback="light"):
    """copy basisset from location LOC of type TYP for ELEMENT w/ no. NUM to FOLDER"""
    rep = f"{num:02d}_{elem}_default"

    try:
        shutil.copy(loc / typ / rep, folder)
    except FileNotFoundError:
        warn(f"{typ} basisset for {elem} not found, use '{fallback}' as fallback")
        shutil.copy(loc / fallback / rep, folder)


def setup_aims(ctx):
    """Set up an aims calculator.

    Parameters
    ----------
    ctx: AimsContext
        The context for the calculation

    Returns
    -------
    calc: ase.calculators.calulator.Calculator
        Calculator object for the calculation
    """

    settings = ctx.settings

    # update k_grid
    if ctx.ref_atoms and "control_kpt" in settings:
        if not "density" in settings.control_kpt:
            warn("'control_kpt' given, but not kpt density. Check!", level=1)
        else:
            kptdensity = settings.control_kpt.density
            k_grid = d2k(ctx.ref_atoms, kptdensity, True)
            talk(f"Update aims k_grid with kpt density of {kptdensity} to {k_grid}")
            ctx.settings.obj["k_grid"] = k_grid
            del ctx.settings["control_kpt"]

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

    if "k_grid" not in aims_settings:
        talk("No k_grid in aims calculator. Check!")

    calc = Aims(**aims_settings)

    return calc
