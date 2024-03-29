"""`vibes run` part of the CLI"""

import click

from vibes.helpers import talk

from .misc import AliasedGroup, complete_files, default_context_settings

# paths = click.Path(exists=True)
paths = complete_files
_prefix = "vibes.run"


@click.command(cls=AliasedGroup)
def run():
    """Run a vibes workflow"""


@run.command()
@click.argument("file", default="aims.in", type=paths)
@click.option("--workdir", help="working directory")
@click.pass_obj
def singlepoint(obj, file, workdir):
    """Run singlepoint calculations from FILE (default: aims.in)"""
    from vibes.calculator import CalculatorContext, run_aims
    from vibes.settings import Settings

    ctx = CalculatorContext(Settings(settings_file=file), workdir=workdir)

    if obj.verbose > 0:
        msg = f"run singlepoint calculations with settings from {file}\n"
        talk(msg, prefix=_prefix)

    run_aims(ctx)


@run.command(context_settings=default_context_settings)
@click.argument("file", default="phonopy.in", type=paths)
@click.option("--workdir", help="work directory")
@click.option("--dry", is_flag=True, help="just prepare inputs in the workdir")
@click.pass_obj
def phonopy(obj, file, workdir, dry):
    """Run a phonopy calculation from FILE (default: phonopy.in)"""
    from vibes.phonopy.context import PhonopyContext
    from vibes.settings import Settings

    ctx = PhonopyContext(Settings(settings_file=file), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run phonopy workflow with settings from {file}\n", prefix=_prefix)

    ctx.run(dry=dry)


@run.command()
@click.argument("file", default="phono3py.in", type=paths)
@click.option("--workdir", help="work directory")
@click.option("--dry", is_flag=True, help="just prepare inputs in the workdir")
@click.pass_obj
def phono3py(obj, file, workdir, dry):
    """Run a phono3py calculation from FILE (default: phono3py.in)"""
    from vibes.phono3py.context import Phono3pyContext
    from vibes.settings import Settings

    ctx = Phono3pyContext(Settings(settings_file=file), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run phono3py workflow with settings from {file}\n", prefix=_prefix)

    ctx.run(dry=dry)


@run.command()
@click.argument("file", default="md.in", type=paths)
@click.option("--workdir", help="working directory")
@click.option("--timeout", default=None, type=int, hidden=True)
@click.pass_obj
def md(obj, file, workdir, timeout):
    """Run an MD simulation from FILE (default: md.in)"""
    from vibes.molecular_dynamics.context import MDContext
    from vibes.settings import Settings

    ctx = MDContext(Settings(settings_file=file), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run MD workflow with settings from {file}\n", prefix=_prefix)

    ctx.run(timeout=timeout)


@run.command()
@click.argument("file", default="relaxation.in", type=paths)
@click.option("--workdir", help="working directory")
@click.option("--timeout", default=None, type=int, hidden=True)
@click.pass_obj
def relaxation(obj, file, workdir, timeout):
    """Run an relaxation from FILE (default: relaxation.in)"""
    from vibes.relaxation.context import RelaxationContext
    from vibes.settings import Settings

    ctx = RelaxationContext(Settings(settings_file=file), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run relaxation workflow with settings from {file}\n", prefix=_prefix)

    ctx.run(timeout=timeout)
