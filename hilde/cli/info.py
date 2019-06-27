"""`hilde info` backend"""

import click

from ase.io import read
from hilde.structure.io import inform
from hilde.scripts.md_sum import md_sum
from hilde.scripts.hilde_phonopy import preprocess

from .misc import AliasedGroup


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
    """inform about a structure in a geometry input file

    Parameters
    ----------
    obj: CliTracker
        The click context passed as an object
    filename: str
        Filename for the geometry file (default: geometry.in)
    format: str
        The format of the geometry file (default: aims)
    symprec: float
        symmetry precision for space group/symmetry recognition
    """

    obj.geometry_file = filename

    atoms = read(obj.geometry_file, format=format)

    inform(atoms, symprec=symprec)


@info.command("settings")
@click.argument("filename", default="settings.in")
@click.pass_obj
def settings_info(obj, filename):
    """inform about content of a settings.in file

    Parameters
    ----------
    filename: str
        The file name for the settings file (default: settings.in)
    """

    obj.settings_file = filename
    click.echo(obj.settings_file.read_text())


@info.command("md")
@click.argument("filename", default="trajectory.son")
@click.option("-p", "--plot", is_flag=True, help="plot a summary")
@click.option("--avg", default=100, help="window size for running avg")
@click.option("-v", "--verbose", is_flag=True, help="be verbose")
def md_info(filename, plot, avg, verbose):
    """inform about content of a settings.in file

    Parameters
    ----------
    filename: str
        Filename of the molecular dynamics trajectory (default: trajectory.son)
    plot: bool
        If True plot a summary
    avg: int
        Window size for the rolling average
    verbose: bool
        If True be verbose
    """

    md_sum(filename, plot, avg, verbose)


@info.command("phonopy")
@click.argument("filename", default="phonopy.in")
def phonopy_info(filename):
    """inform about a phonopy calculation before it is started

    Parameters
    ----------
    filename: str
        Filename for the phonopy settings file (default: phonopy.in)
    """

    preprocess(filename=None, settings_file=filename, dimension=None, format=None)
