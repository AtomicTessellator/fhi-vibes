"""`hilde info` backend"""

from pathlib import Path

import click

from ase.io import read
from hilde.structure.io import inform
from hilde.scripts.md_sum import md_sum
from hilde.scripts.hilde_phonopy import preprocess

from .misc import AliasedGroup


def check_file(filename):
    """check if file exists"""
    if not Path(filename).exists():
        raise click.FileError(filename, hint="not found")


@click.command(cls=AliasedGroup)
def info():
    """inform about content of a file"""
    pass


@info.command("geometry")
@click.argument("filename", default="geometry.in")
@click.option("--format", default="aims", show_default=True)
@click.option("-t", "--symprec", default=1e-5, show_default=True)
@click.pass_obj
def geometry_info(obj, filename, format, symprec):
    """inform about a structure in a geometry input file"""

    check_file(filename)

    atoms = read(filename, format=format)

    inform(atoms, symprec=symprec)


@info.command("md")
@click.argument("filename", default="trajectory.son")
@click.option("-p", "--plot", is_flag=True, help="plot a summary")
@click.option("--avg", default=100, help="window size for running avg")
@click.option("-v", "--verbose", is_flag=True, help="be verbose")
def md_info(filename, plot, avg, verbose):
    """inform about content of a settings.in file"""

    check_file(filename)

    md_sum(filename, plot, avg, verbose)


@info.command("phonopy")
@click.argument("filename", default="phonopy.in")
def phonopy_info(filename):
    """inform about a phonopy calculation before it is started"""

    check_file(filename)
    preprocess(filename=None, settings_file=filename, dimension=None, format=None)
