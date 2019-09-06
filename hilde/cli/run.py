"""`hilde run` part of the CLI"""

import click

from hilde.aims.context import AimsContext
from hilde.phonopy.context import PhonopyContext
from hilde.molecular_dynamics.context import MDContext
from hilde.settings import Settings
from hilde.helpers import talk

from .misc import AliasedGroup, complete_filenames

# paths = click.Path(exists=True)
paths = complete_filenames
_prefix = "hilde.run"


@click.command(cls=AliasedGroup)
def run():
    """run a hilde workflow"""


@run.command("aims")
@click.argument("settings", default="aims.in", type=paths)
@click.option("--workdir", help="working directory")
@click.pass_obj
def aims_run(obj, settings, workdir):
    """run one or several aims calculations from SETTINGS (default: aims.in)"""
    from hilde.aims.workflow import run_aims

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
    from hilde.phonopy.workflow import run_phonopy

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
    """run an MD simulation from SETTINS (default: md.in)"""
    ctx = MDContext(Settings(settings_file=settings), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run MD workflow with settings from {settings}\n", prefix=_prefix)

    ctx.run(timeout=timeout)
