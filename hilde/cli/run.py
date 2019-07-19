"""`hilde run` part of the CLI"""

import click

from hilde.aims.context import AimsContext
from hilde.phonopy.context import PhonopyContext
from hilde.molecular_dynamics.context import MDContext
from hilde.settings import Settings
from hilde.helpers import cwd, talk

from .misc import AliasedGroup, complete_filenames

# paths = click.Path(exists=True)
paths = complete_filenames
_prefix = "hilde.run"


@click.command(cls=AliasedGroup)
def run():
    """run a hilde workflow"""


@run.command("aims")
@click.option("--workdir", help="working directory")
@click.option("--settings", default="aims.in", show_default=True, type=paths)
@click.pass_obj
def aims_run(obj, workdir, settings):
    """run one or several aims calculations"""
    from hilde.aims.workflow import run_aims

    ctx = AimsContext(Settings(settings_file=settings), workdir=workdir)

    run_aims(ctx)


@run.command("phonopy")
@click.option("--workdir", help="work directory")
@click.option("--settings", default="phonopy.in", show_default=True, type=paths)
@click.option("--dry", is_flag=True, help="just prepare inputs in the workdir")
@click.pass_obj
def phonopy_run(obj, workdir, settings, dry):
    """run a phonopy calculation"""
    from hilde.phonopy.workflow import run_phonopy

    ctx = PhonopyContext(Settings(settings_file=settings), workdir=workdir)

    if dry:
        with cwd(ctx.workdir, mkdir=True):
            ctx.ref_atoms.write(ctx.settings.geometry.file)
            ctx.settings.obj["workdir"] = "."
            ctx.settings.write()
    else:
        run_phonopy(ctx=ctx)


@run.command("md")
@click.option("--workdir", help="working directory")
@click.option("--settings", default="md.in", show_default=True, type=paths)
@click.option("--timeout", default=None, type=int, hidden=True)
@click.pass_obj
def md_run(obj, workdir, settings, timeout):
    """run an MD simulation"""
    ctx = MDContext(Settings(settings_file=settings), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run MD workflow with settings from {settings}\n", prefix=_prefix)

    ctx.run(timeout=timeout)
