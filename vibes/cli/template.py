"""`vibes input` part of the CLI"""

import click

from vibes import keys
from vibes.templates import config_files, settings

from .misc import AliasedGroup, ClickAliasedGroup

try:
    import importlib.resources as pkg_resources
except ModuleNotFoundError:
    import importlib_resources as pkg_resources


@click.command(cls=AliasedGroup)
@click.option("--allow_overwrite", is_flag=True, show_default=True)
@click.pass_obj
def template(obj, allow_overwrite):
    """provide template input files for tasks and workflows"""
    obj.allow_overwrite = allow_overwrite


@template.command(cls=ClickAliasedGroup)
def calculator():
    """Calculator templates: aims, lj"""


@calculator.command()
@click.argument("file", default="aims.in")
@click.pass_obj
def aims(obj, file):
    """provide template aims.in for aims calculator"""

    print_input(obj, "aims", file)


@calculator.command()
@click.argument("file", default="lj.in")
@click.pass_obj
def lj(obj, file):
    """provide template lj.in for Lennard-Jones calculator for solid Argon"""

    print_input(obj, "lj", file)


@template.command()
@click.argument("file", default="phonopy.in")
@click.pass_obj
def phonopy(obj, file):
    """provide template phonopy.in for phonopy workflow."""
    from vibes.phonopy.context import PhonopyContext

    ctx = PhonopyContext()
    ctx.settings.print()


@template.command()
@click.argument("file", default="md.in")
@click.option("--nvt", is_flag=True)
@click.option("--npt", is_flag=True)
@click.pass_obj
def md(obj, file, nvt, npt):
    """provide template md.in for molecular dynamics workflow."""
    from vibes.molecular_dynamics.context import MDContext

    if nvt:
        ensemble = keys.nvt
    elif npt:
        ensemble = keys.npt
    else:
        ensemble = keys.nve

    ctx = MDContext(ensemble=ensemble)
    ctx.settings.print()


@template.command()
@click.argument("file", default="relaxation.in")
@click.pass_obj
def relaxation(obj, file):
    """provide template relaxation.in for relaxation workflow."""
    from vibes.relaxation.context import RelaxationContext

    ctx = RelaxationContext()
    ctx.settings.print()


@template.command(cls=ClickAliasedGroup)
def configuration():
    """Configuration templates: .vibesrc, .fireworksrc"""


@configuration.command()
@click.argument("file", default="vibesrc")
@click.pass_obj
def vibes(obj, file):
    """provide template .vibesrc for the configuration"""

    print_input(obj, "vibesrc.template", file, from_folder=config_files)


@configuration.command()
@click.argument("file", default="fireworksrc")
@click.pass_obj
def fireworks(obj, file):
    """provide template fireworksrc.template for the configuration"""

    print_input(obj, "fireworksrc.template", file, from_folder=config_files)


@template.command()
@click.argument("file", default="slurm.in")
@click.pass_obj
def slurm(obj, file):
    """provide template slurm settings"""

    print_input(obj, "slurm.in", file, from_folder=config_files)


def print_input(obj, name, file, from_folder=settings):
    """write the input function"""

    input_file = pkg_resources.read_text(from_folder, name)

    click.echo(input_file)
