"""`vibes run` part of the CLI"""

import click

from vibes.helpers import talk

from .misc import AliasedGroup, complete_files

# paths = click.Path(exists=True)
paths = complete_files
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
    from vibes.calculator import CalculatorContext, run_aims

    ctx = CalculatorContext(Settings(settings_file=settings), workdir=workdir)

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

    ctx = PhonopyContext(Settings(settings_file=settings), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run phonopy workflow with settings from {settings}\n", prefix=_prefix)

    ctx.run(dry=dry)


@run.command("phono3py")
@click.argument("settings", default="phono3py.in", type=paths)
@click.option("--workdir", help="work directory")
@click.option("--kappa", help="")
@click.option("--dry", is_flag=True, help="just prepare inputs in the workdir")
@click.pass_obj
def phono3py_run(obj, settings, workdir, dry):
    """run a phonopy calculation for SETTINGS (default: phonopy.in)"""
    from vibes.settings import Settings
    from vibes.phono3py.context import Phono3pyContext

    ctx = Phono3pyContext(Settings(settings_file=settings), workdir=workdir)

    if obj.verbose > 0:
        talk(f"run phono3py workflow with settings from {settings}\n", prefix=_prefix)

    ctx.run(dry=dry)


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
