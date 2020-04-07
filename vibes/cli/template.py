"""`vibes input` part of the CLI"""

from pathlib import Path

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

    write_input(obj, "aims", file)


@template.command("phonopy")
@click.argument("file", default="phonopy.in")
@click.pass_obj
def phonopy_input(obj, file):
    """provide template phonopy.in for phonopy workflow."""
    from vibes.phonopy.context import PhonopyContext

    ctx = PhonopyContext()
    write_settings(obj, ctx.settings, file)


@template.command("md")
@click.argument("file", default="md.in")
@click.pass_obj
def md_input(obj, file):
    """provide template md.in for molecular dynamics workflow."""
    from vibes.molecular_dynamics.context import MDContext

    ctx = MDContext()
    write_settings(obj, ctx.settings, file)


@template.command("relaxation")
@click.argument("file", default="relaxation.in")
@click.pass_obj
def relaxation_input(obj, file):
    """provide template relaxation.in for relaxation workflow."""
    from vibes.relaxation.context import RelaxationContext

    ctx = RelaxationContext()
    write_settings(obj, ctx.settings, file)


@template.command("configuration")
@click.argument("file", default="vibesrc")
@click.pass_obj
def configuration_input(obj, file):
    """provide template vibesrc.template for the configuration"""

    write_input(obj, "vibesrc.template", file, from_folder=config_files)


@template.command("fireworks_configuration")
@click.argument("file", default="fireworksrc")
@click.pass_obj
def fireworks_configuration_input(obj, file):
    """provide template fireworksrc.template for the configuration"""

    write_input(obj, "fireworksrc.template", file, from_folder=config_files)


@template.command("slurm")
@click.argument("file", default="slurm.in")
@click.pass_obj
def slurm_input(obj, file):
    """provide template slurm settings"""

    write_input(obj, "slurm.in", file, from_folder=config_files)


def write_input(obj, name, file, from_folder=settings):
    """write the input function"""

    input_file = pkg_resources.read_text(from_folder, name)

    outfile = Path(file)

    if not obj.allow_overwrite and outfile.exists():
        msg = f"{outfile} exists."
        raise click.ClickException(msg)

    outfile.write_text(input_file)

    click.echo(f"Default {name} settings file written to {file}.")


def write_settings(obj, settings, file):
    """write the settings"""

    outfile = Path(file)

    if outfile.exists() and obj.allow_overwrite:
        click.echo(f"Remove `{outfile}`")
        outfile.unlink()

    settings.write(outfile)

    click.echo(f"Settings written to `{outfile}`")
