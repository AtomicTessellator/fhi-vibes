"""`hilde run` part of the CLI"""

import click

from hilde.aims.context import AimsContext
from hilde.phonopy.context import PhonopyContext
from hilde.molecular_dynamics.context import MDContext
from hilde.settings import Settings

from .misc import AliasedGroup, click, complete_filenames

# paths = click.Path(exists=True)
paths = complete_filenames


@click.command(cls=AliasedGroup)
def run():
    """run a hilde workflow"""


@run.command("aims")
@click.option("--workdir", default="aims", help="working directory")
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
@click.pass_obj
def phonopy_run(obj, workdir, settings):
    """run a phonopy calculation"""
    from hilde.phonopy.workflow import run_phonopy

    ctx = PhonopyContext(Settings(settings_file=settings), workdir=workdir)

    run_phonopy(ctx=ctx)


@run.command("md")
@click.option("--workdir", default="md", help="working directory")
@click.option("--settings", default="md.in", show_default=True, type=paths)
@click.pass_obj
def md_run(obj, workdir, settings):
    """run an MD simulation"""
    from hilde.molecular_dynamics.workflow import run_md

    ctx = MDContext(Settings(settings_file=settings), workdir=workdir)

    if obj.verbose > 0:
        click.echo(f"run MD workflow with settings from {settings}\n")
    run_md(ctx=ctx)
