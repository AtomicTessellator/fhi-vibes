"""`hilde input` part of the CLI"""

from pathlib import Path
import importlib.resources as pkg_resources

import click
from hilde.templates import settings, config_files

from .misc import AliasedGroup


@click.command(cls=AliasedGroup)
@click.option("--full", is_flag=True, help="list more options", show_default=True)
@click.option("--allow_overwrite", is_flag=True, show_default=True)
@click.pass_obj
def template(obj, full, allow_overwrite):
    """provide template input files for tasks and workflows"""
    obj.full_input = full
    obj.allow_overwrite = allow_overwrite


@template.command("modify")
@click.argument("filename", default="settings.in")
@click.pass_obj
def modify_input(obj, filename):
    """modify an input file"""

    click.echo("please come back later")


@template.command("aims")
@click.argument("filename", default="aims.in")
@click.pass_obj
def aims_input(obj, filename):
    """provide template settings.in for aims calculation"""

    write_input(obj, "aims", filename)


@template.command("phonopy")
@click.argument("filename", default="phonopy.in")
@click.pass_obj
def phonopy_input(obj, filename):
    """provide template phonopy.in for phonopy workflow."""

    write_input(obj, "phonopy", filename)


@template.command("md")
@click.argument("filename", default="md.in")
@click.pass_obj
def md_input(obj, filename):
    """provide template md.in for molecular dynamics workflow."""

    write_input(obj, "md", filename)


@template.command("configuration")
@click.argument("filename", default="hilderc")
@click.pass_obj
def configuration_input(obj, filename):
    """provide template hilderc.template for the configuration"""

    write_input(obj, "hilderc.template", filename, from_folder=config_files)


def write_input(obj, name, filename, from_folder=settings):
    """write the input function"""

    if obj.full_input:
        name += "_full"

    input_file = pkg_resources.read_text(from_folder, name)

    outfile = Path(filename)

    if not obj.allow_overwrite and outfile.exists():
        msg = f"{outfile} exists."
        raise click.ClickException(msg)

    outfile.write_text(input_file)

    click.echo(f"Default {name} settings file written to {filename}.")
