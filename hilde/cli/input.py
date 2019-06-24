"""`hilde input` part of the CLI"""

from pathlib import Path
import importlib.resources as pkg_resources

import click
from hilde.templates import settings

from .misc import AliasedGroup


@click.command(cls=AliasedGroup)
@click.option("--full", is_flag=True, help="list more options", show_default=True)
@click.option("--allow_overwrite", is_flag=True, show_default=True)
@click.pass_obj
def input(obj, full, allow_overwrite):
    """provide template input files for tasks and workflows"""
    obj.full_input = full
    obj.allow_overwrite = allow_overwrite


@input.command("modify")
@click.argument("filename", default="settings.in")
@click.pass_obj
def modify_input(obj, filename):
    """modify an input file"""

    click.echo("please come back later")


@input.command("aims")
@click.argument("filename", default="aims.in")
@click.pass_obj
def aims_input(obj, filename):
    """provide template settings.in for aims calculation"""

    write_input("aims", filename, obj.allow_overwrite)


@input.command("phonopy")
@click.argument("filename", default="phonopy.in")
@click.pass_obj
def phonopy_input(obj, filename):
    """provide template phonopy.in for phonopy workflow."""

    write_input("phonopy", filename, obj.allow_overwrite)


@input.command("md")
@click.argument("filename", default="md.in")
@click.pass_obj
def md_input(obj, filename):
    """provide template md.in for molecular dynamics workflow."""

    write_input("md", filename, obj.allow_overwrite)


@click.pass_obj
def write_input(obj, name, filename, allow_overwrite):
    """write the input function"""

    if obj.full_input:
        name += "_full"

    template = pkg_resources.read_text(settings, name)

    outfile = Path(filename)

    if not allow_overwrite and outfile.exists():
        msg = f"{outfile} exists."
        raise click.ClickException(msg)

    outfile.write_text(template)

    click.echo(f"Default {name} settings file written to {filename}.")
