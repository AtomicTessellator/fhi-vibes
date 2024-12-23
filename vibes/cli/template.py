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
def template():
    """Provide template input files for tasks and workflows"""


@template.command(cls=ClickAliasedGroup)
def calculator():
    """Calculator templates: aims, lj"""


@calculator.command()
def aims():
    """Provide template input for aims calculator"""
    print_input("aims")


@calculator.command()
def lj():
    """Provide template input for Lennard-Jones calculator for solid Argon"""
    print_input("lj")


@template.command()
def phonopy():
    """Provide template input for phonopy workflow."""
    from vibes.phonopy.context import PhonopyContext

    ctx = PhonopyContext()
    ctx.settings.print()

@template.command()
def phono3py():
    """Provide template input for phono3py workflow."""
    from vibes.phono3py.context import Phono3pyContext

    ctx = Phono3pyContext()
    ctx.settings.print()

@template.command()
@click.option("--nvt", is_flag=True, help="Use Langevin thermostat for NVT simulation")
@click.option("--npt", is_flag=True, help="Use Berendsen algorithm for NPT simulation")
def md(nvt, npt):
    """Provide template input for MD simulation (default: NVE)"""
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
def relaxation():
    """Provide template input for relaxation workflow."""
    from vibes.relaxation.context import RelaxationContext

    ctx = RelaxationContext()
    ctx.settings.print()


@template.command(cls=ClickAliasedGroup)
def configuration():
    """Configuration templates: .vibesrc, .fireworksrc"""


@configuration.command()
def vibes():
    """Provide template input for .vibesrc"""
    print_input("vibesrc.template", from_folder=config_files)


@configuration.command()
def fireworks():
    """Provide template inpurt for .fireworksrc"""
    print_input("fireworksrc.template", from_folder=config_files)


@template.command()
def slurm():
    """Provide template slurm settings"""
    print_input("slurm.in", from_folder=config_files)


def print_input(name, from_folder=settings):
    """Write the input function"""
    input_file = pkg_resources.read_text(from_folder, name)

    click.echo(input_file)
