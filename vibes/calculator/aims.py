"""FHI-aims related setup"""
import shutil
from pathlib import Path

# from ase.calculators.aims import Aims
from ase.calculators.aims import Aims

from vibes.helpers.k_grid import d2k
from vibes.helpers.warnings import warn

from ._defaults import basisset_choices, basisset_default, basisset_key, name, talk
from .context import CalculatorContext

_fallback = "light"


class BasissetError(RuntimeError):
    """Raise when the basisset was set up incorrectly"""


def create_species_dir(
    ctx: CalculatorContext, folder: str = "basissets", fallback: str = _fallback
) -> str:
    """ create a custom bassiset folder for the computation

    Args:
        ctx: The context for the calculation
        folder: Folder to store the basisset

    Returns:
        The absolute file path to the species directory

    """
    loc = ctx.basisset_location
    settings = ctx.settings

    # if old section with `basisset.type` is used:
    if basisset_key not in settings:
        msg = "basissets not specified in settings.file."
        raise BasissetError(msg)

    default = basisset_default
    fallback = basisset_default

    if "default" in settings.basissets:
        default = settings.basissets.default
    if "fallback" in settings.basissets:
        fallback = settings.basissets.fallback

    for key in (default, fallback):
        if key not in basisset_choices:
            raise BasissetError(f"Species default '{key}' unknown.")

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

    for (key, val) in key_vals:
        # copy the respective basisset
        add_basisset(loc, val, key, dct[key], folder, fallback=fallback)
        del dct[key]

    # add remaining ones
    for key in dct.keys():
        # copy the respective basisset
        add_basisset(loc, default, key, dct[key], folder, fallback=fallback)

    return str(folder.absolute())


def add_basisset(
    loc: str,
    typ: str,
    elem: str,
    num: int,
    folder: str,
    fallback: str = _fallback,
    verbose: bool = True,
):
    """copy basisset from location LOC of type TYP for ELEMENT w/ no. NUM to FOLDER"""
    rep = f"{num:02d}_{elem}_default"

    msg = f"Add basisset `{typ}` for atom `{elem}` to basissets folder."
    talk(msg, verbose=verbose)

    try:
        shutil.copy(loc / typ / rep, folder)
    except FileNotFoundError:
        warn(f"{typ} basisset for {elem} not found, use '{fallback}' as fallback")
        shutil.copy(loc / fallback / rep, folder)


def setup_aims(
    ctx: CalculatorContext, verbose: bool = True, make_species_dir: bool = True
) -> Aims:
    """Set up an aims calculator.

    Args:
        ctx (AimsContext): The context for the calculation
        verbose (bool): inform about the calculator details

    Returns:
        Calculator object for the calculation

    """
    settings = ctx.settings

    # update k_grid
    if ctx.ref_atoms and "control_kpt" in settings:
        if "density" not in settings.control_kpt:
            warn("'control_kpt' given, but not kpt density. Check!", level=1)
        else:
            kptdensity = settings.control_kpt.density
            k_grid = d2k(ctx.ref_atoms, kptdensity, True)
            talk(f"Update aims k_grid with kpt density of {kptdensity} to {k_grid}")
            ctx.settings.obj["k_grid"] = k_grid
            del ctx.settings["control_kpt"]

    aims_settings = settings.obj

    ase_settings = {"aims_command": settings.machine.aims_command}

    if "socketio" in settings:
        host = settings.socketio.get("host", "localhost")
        port = settings.socketio.port

        if settings.socketio.get("unixsocket", None) is not None:
            host = f"UNIX:{settings.socketio.unixsocket}"
            port = settings.socketio.get("port", 31415)

        if port is not None:
            aims_settings.update({"use_pimd_wrapper": (host, port)})

    # create basissetfolder
    if make_species_dir:
        species_dir = create_species_dir(ctx)
    else:
        species_dir = str(ctx.basisset_location / settings.basissets.default)

    ase_settings["species_dir"] = species_dir

    aims_settings = {**aims_settings, **ase_settings}

    if "k_grid" not in aims_settings:
        talk("No k_grid in aims calculator. Check!")

    talk(f"Calculator: {name}", verbose=verbose)
    msg = ["settings:", *[f"  {k}: {v}" for k, v in aims_settings.items()]]
    talk(msg, verbose=verbose)

    calculator = Aims(**aims_settings)

    return calculator
