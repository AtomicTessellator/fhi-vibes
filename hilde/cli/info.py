"""`hilde info` backend"""

from ase.io import read
from hilde.structure.io import inform
from hilde.scripts.md_sum import md_sum
from hilde.scripts.hilde_phonopy import preprocess

from .misc import AliasedGroup, click, complete_filenames


@click.command(cls=AliasedGroup)
def info():
    """inform about content of a file"""


@info.command("geometry")
@click.argument("filename", default="geometry.in", type=complete_filenames)
@click.option("--format", default="aims", show_default=True)
@click.option("-t", "--symprec", default=1e-5, show_default=True)
@click.option("-v", "--verbose", is_flag=True)
@click.pass_obj
def geometry_info(obj, filename, format, symprec, verbose):
    """inform about a structure in a geometry input file"""

    atoms = read(filename, format=format)

    verbosity = 1
    if verbose:
        verbosity = 2

    inform(atoms, symprec=symprec, verbosity=verbosity)


@info.command("settings")
@click.argument("filename", default="settings.in")
@click.pass_obj
def settings_info(obj, filename):
    """inform about content of a settings.in file"""

    obj.settings_file = filename
    click.echo(obj.settings_file.read_text())


@info.command("md")
@click.argument("filename", default="trajectory.son", type=complete_filenames)
@click.option("-p", "--plot", is_flag=True, help="plot a summary")
@click.option("--avg", default=100, help="window size for running avg")
@click.option("-v", "--verbose", is_flag=True, help="be verbose")
def md_info(filename, plot, avg, verbose):
    """inform about content of a settings.in file"""

    md_sum(filename, plot, avg, verbose)


@info.command("phonopy")
@click.argument("filename", default="phonopy.in", type=complete_filenames)
def phonopy_info(filename):
    """inform about a phonopy calculation before it is started"""

    preprocess(filename=None, settings_file=filename, dimension=None, format=None)
