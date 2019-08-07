""" Provide a full highlevel phonopy workflow

    Input: geometry.in and settings.in
    Output: geometry.in.supercell and trajectory.son """

from hilde.tasks import calculate_socket
from hilde.helpers.restarts import restart

from hilde.aims.context import AimsContext
from hilde.aims.setup import setup_aims
from hilde.helpers import talk

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


def bootstrap(ctx):
    """load settings, prepare atoms, calculator, and phonopy

    Args:
        ctx (PhonopyContext): The context for the calculation

    Returns:
        dict: The necessary information to run the workflow with the following items
            atoms_to_calculate (list): list of the displaced supercells
            calculator (ase.calculators.calulator.Calculator):use to calculate forces
            metadata (dict): metadata for the phonon calculation
            workdir (str or Path): working directory for the calculation
            settings (Settings): settings for the workflow
    """
    if ctx.name.lower() == "phonopy":
        from hilde.phonopy.wrapper import preprocess
    elif ctx.name.lower() == "phono3py":
        from hilde.phono3py.wrapper import preprocess

    # Phonopy preprocess
    phonon, supercell, scs = preprocess(atoms=ctx.ref_atoms, **ctx.settings.obj)

    # if calculator not given, create an aims context for this calculation
    if ctx.settings.atoms and  ctx.settings.atoms.calc:
        calc = ctx.settings.atoms.calc
    else:
        aims_ctx = AimsContext(settings=ctx.settings, workdir=ctx.workdir)
        # set reference structure for aims calculation and make sure forces are computed
        aims_ctx.ref_atoms = supercell
        aims_ctx.settings.obj["compute_forces"] = True

        calc = setup_aims(aims_ctx)

    # save metadata
    metadata = metadata2dict(phonon, calc)

    return {
        "atoms_to_calculate": scs,
        "calculator": calc,
        "metadata": metadata,
        "workdir": ctx.workdir,
        "settings": ctx.settings,
        "save_input": True,
        "backup_after_calculation": False,
        **ctx.settings.obj,
    }
