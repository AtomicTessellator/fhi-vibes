""" Provide a full highlevel phonopy workflow

    Input: geometry.in and settings.in
    Output: geometry.in.supercell and trajectory.son """

from hilde.tasks import calculate_socket
from hilde.helpers.restarts import restart

from hilde.aims.context import AimsContext
from hilde.aims.setup import setup_aims
from hilde.helpers import talk

from .context import PhonopyContext
from .postprocess import postprocess
from . import metadata2dict


def run_phonopy(**kwargs):
    """ high level function to run phonopy workflow """

    args = bootstrap(**kwargs)
    workdir = args["workdir"]

    talk(f"Run phonopy workflow in working directory\n  {workdir}")

    try:
        postprocess(**args)
        msg = "** Postprocess could be performed from previous calculations. Check"
        msg += f"\n**  {workdir}"
        exit(msg)
    except (FileNotFoundError, RuntimeError):
        completed = calculate_socket(**args)

    if not completed:
        restart()
    else:
        talk("Start postprocess.")
        postprocess(**args)
        talk("done.")


def bootstrap(ctx=None, name=None, settings=None, workdir=None, **kwargs):
    """load settings, prepare atoms, calculator, and phonopy

    Parameters
    ----------
    ctx: PhonopyContext
        The context for the calculation
    name: str
        Name of the type of calculation
    settings: Settings
        settings for the workflow
    workdir: str or Path
        The working directory for the calculation

    Returns
    -------
    dict
        The necessary information to run the workflow with the following items

        atoms_to_calculate: list of ase.atoms.Atoms
            The list of the displaced supercells
        calculator: ase.calculators.calulator.Calculator
            The calculator used to calculate for forces in each supercell
        metadata: dict
            The metadata for the phonon calculation
        workdir: str or Path
            The working directory for the calculation
        settings: Settings
            The settings for the workflow

        Additional key/value pairs in settings.obj
    """
    if ctx is None:
        ctx = PhonopyContext(settings=settings)
    if workdir:
        ctx.workdir = workdir

    settings = ctx.settings

    if not name:
        name = ctx.name

    if name.lower() == "phonopy":
        from hilde.phonopy.wrapper import preprocess
    elif name.lower() == "phono3py":
        from hilde.phono3py.wrapper import preprocess

    # Phonopy preprocess
    phonon, supercell, scs = preprocess(atoms=ctx.ref_atoms, **ctx.settings.obj)

    # if calculator not given, create an aims context for this calculation
    if "calculator" not in kwargs:
        aims_ctx = AimsContext(settings=ctx.settings, workdir=ctx.workdir)
        # set reference structure for aims calculation and make sure forces are computed
        aims_ctx.ref_atoms = supercell
        aims_ctx.settings.obj["compute_forces"] = True

        calc = setup_aims(aims_ctx)
    else:
        calc = kwargs["calculator"]

    # save metadata
    metadata = metadata2dict(phonon, calc)

    return {
        "atoms_to_calculate": scs,
        "calculator": calc,
        "metadata": metadata,
        "workdir": ctx.workdir,
        "settings": settings,
        "backup_after_calculation": False,
        **settings.obj,
    }
