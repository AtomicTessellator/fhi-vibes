"""`hilde run` part of the CLI"""

import click

from hilde.aims.context import AimsContext
from hilde.phonopy.context import PhonopyContext
from hilde.molecular_dynamics.context import MDContext

from .misc import AliasedGroup


@click.command(cls=AliasedGroup)
def run():
    """run a hilde workflow"""
    pass


@run.command("aims")
@click.option("--workdir", default="aims", help="working directory")
@click.option("--settings", default="aims.in", show_default=True)
@click.pass_obj
def aims_run(obj, workdir, settings):
    """run and aims calculation"""
    from hilde.aims.workflow import run_aims

    context = AimsContext(settings_file=settings, workdir=workdir)

    run_aims(context)


@run.command("phonopy")
@click.option("--workdir", default="phonopy", show_default=True, help="work directory")
@click.option("--settings", default="phonopy.in", show_default=True)
@click.pass_obj
def phonopy_run(obj, workdir, settings):
    """run and aims calculation"""
    from hilde.phonopy.workflow import run_phonopy

    ctx = PhonopyContext(settings_file=settings, workdir=workdir)

    run_phonopy(ctx=ctx)


@run.command("md")
@click.option("--workdir", default="md", help="working directory")
@click.option("--settings", default="md.in", show_default=True)
@click.pass_obj
def md_run(obj, workdir, settings):
    """run and aims calculation"""
    from hilde.molecular_dynamics.workflow import run_md

    ctx = MDContext(settings_file=settings, workdir=workdir)

    if obj.verbose > 0:
        click.echo(f"run MD workflow with settings from {settings}\n")
    run_md(ctx=ctx)
