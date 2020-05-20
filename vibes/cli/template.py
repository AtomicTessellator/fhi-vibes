"""`vibes input` part of the CLI"""

import click

from vibes.templates import config_files, settings

from .misc import AliasedGroup

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


@template.command("aims")
@click.argument("file", default="aims.in")
@click.pass_obj
def aims_input(obj, file):
    """provide template settings.in for aims calculation"""

    print_input(obj, "aims", file)


@template.command("phonopy")
@click.argument("file", default="phonopy.in")
@click.pass_obj
def phonopy_input(obj, file):
    """provide template phonopy.in for phonopy workflow."""
    from vibes.phonopy.context import PhonopyContext

    ctx = PhonopyContext()
    ctx.settings.print()


@template.command("md")
@click.argument("file", default="md.in")
@click.pass_obj
def md_input(obj, file):
    """provide template md.in for molecular dynamics workflow."""
    from vibes.molecular_dynamics.context import MDContext

    ctx = MDContext()
    ctx.settings.print()


@template.command("relaxation")
@click.argument("file", default="relaxation.in")
@click.pass_obj
def relaxation_input(obj, file):
    """provide template relaxation.in for relaxation workflow."""
    from vibes.relaxation.context import RelaxationContext

    ctx = RelaxationContext()
    ctx.settings.print()


@template.command("configuration")
@click.argument("file", default="vibesrc")
@click.pass_obj
def configuration_input(obj, file):
    """provide template vibesrc.template for the configuration"""

    print_input(obj, "vibesrc.template", file, from_folder=config_files)


@template.command("fireworks_configuration")
@click.argument("file", default="fireworksrc")
@click.pass_obj
def fireworks_configuration_input(obj, file):
    """provide template fireworksrc.template for the configuration"""

    print_input(obj, "fireworksrc.template", file, from_folder=config_files)


@template.command("slurm")
@click.argument("file", default="slurm.in")
@click.pass_obj
def slurm_input(obj, file):
    """provide template slurm settings"""

    print_input(obj, "slurm.in", file, from_folder=config_files)


def print_input(obj, name, file, from_folder=settings):
    """write the input function"""

    input_file = pkg_resources.read_text(from_folder, name)

    click.echo(input_file)
