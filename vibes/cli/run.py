"""`vibes run` part of the CLI"""

import click

# from vibes.aims.context import AimsContext
# from vibes.phonopy.context import PhonopyContext
# from vibes.molecular_dynamics.context import MDContext
# from vibes.settings import Settings
from vibes.helpers import talk

from .misc import AliasedGroup, complete_filenames


# paths = click.Path(exists=True)
paths = complete_filenames
_prefix = "vibes.run"


@click.command(cls=AliasedGroup)
def run():
    """run a vibes workflow"""


@run.command("aims")
@click.argument("settings", default="aims.in", type=paths)
@click.option("--workdir", help="working directory")
@click.pass_obj
def aims_run(obj, settings, workdir):
    """run one or several aims calculations from SETTINGS (default: aims.in)"""
    from vibes.settings import Settings
    from vibes.aims import AimsContext, run_aims

    ctx = AimsContext(Settings(settings_file=settings), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run aims calculations with settings from {settings}\n", prefix=_prefix)

    run_aims(ctx)


@run.command("phonopy")
@click.argument("settings", default="phonopy.in", type=paths)
@click.option("--workdir", help="work directory")
@click.option("--dry", is_flag=True, help="just prepare inputs in the workdir")
@click.pass_obj
def phonopy_run(obj, settings, workdir, dry):
    """run a phonopy calculation for SETTINGS (default: phonopy.in)"""
    from vibes.settings import Settings
    from vibes.phonopy.context import PhonopyContext
    from vibes.phonopy.workflow import run_phonopy

    ctx = PhonopyContext(Settings(settings_file=settings), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run phonpoy workflow with settings from {settings}\n", prefix=_prefix)

    run_phonopy(ctx=ctx, dry=dry)


@run.command("md")
@click.argument("settings", default="md.in", type=paths)
@click.option("--workdir", help="working directory")
@click.option("--timeout", default=None, type=int, hidden=True)
@click.pass_obj
def md_run(obj, settings, workdir, timeout):
    """run an MD simulation from SETTINGS (default: md.in)"""
    from vibes.settings import Settings
    from vibes.molecular_dynamics.context import MDContext

    ctx = MDContext(Settings(settings_file=settings), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run MD workflow with settings from {settings}\n", prefix=_prefix)

    ctx.run(timeout=timeout)


@run.command("relaxation")
@click.argument("settings", default="relaxation.in", type=paths)
@click.option("--workdir", help="working directory")
@click.option("--timeout", default=None, type=int, hidden=True)
@click.pass_obj
def relaxation_run(obj, settings, workdir, timeout):
    """run an relaxation from SETTINGS (default: relaxation.in)"""
    from vibes.settings import Settings
    from vibes.relaxation.context import RelaxationContext

    ctx = RelaxationContext(Settings(settings_file=settings), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run relaxation workflow with settings from {settings}\n", prefix=_prefix)

    ctx.run(timeout=timeout)