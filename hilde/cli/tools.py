"""hilde CLI tools"""

import click

from hilde.scripts.get_relaxation_info import get_relaxation_info
from hilde.scripts.refine_geometry import refine_geometry
from hilde.scripts.make_supercell import make_supercell
from hilde.scripts.create_samples import create_samples
from hilde.scripts.suggest_k_grid import suggest_k_grid
from hilde.scripts.remap_phonopy_forceconstants import remap_phonopy_force_constants
from hilde.scripts.nomad_upload import nomad_upload

from .misc import AliasedGroup


@click.command(cls=AliasedGroup)
def tools():
    """tools"""


@tools.command(cls=AliasedGroup)
def geometry():
    """tools for working with structures"""
    ...


@geometry.command("refine")
@click.argument("filename", default="geometry.in")
@click.option("-prim", "--primitive", is_flag=True)
@click.option("-conv", "--conventional", is_flag=True)
@click.option("--center", is_flag=True)
@click.option("--origin", is_flag=True)
@click.option("-cart", "--cartesian", is_flag=True)
@click.option("--format", default="aims")
@click.option("-t", "--symprec", default=1e-5)
def geometry_refine(*args, **kwargs):
    """hilde.scripts.refine_geometry"""
    refine_geometry(*args, **kwargs)


@tools.command("make_supercell")
@click.argument("filename", default="geometry.in")
@click.option("-d", "--dimension", type=float)
@click.option("-n", "--n_target", type=int)
@click.option("--deviation", default=0.2)
@click.option("--dry", is_flag=True)
@click.option("--format", default="aims")
@click.option("--scaled", is_flag=True)
def tool_make_supercell(filename, dimension, n_target, deviation, dry, format, scaled):
    """create a supercell of desired shape or size"""
    make_supercell(filename, dimension, n_target, deviation, dry, format, scaled)


# @click.command(cls=AliasedGroup)
# def utils():
#     """utilities, for example `aims get_relaxation_info`"""


@tools.command("get_relaxation_info")
@click.argument("filenames", nargs=-1)
def relaxation_info(filenames):
    """analyize aims relaxation"""
    get_relaxation_info(filenames)


@tools.command("create_samples")
@click.argument("filename")
@click.option("-T", "--temperature", type=float, help="Temperature in Kelvin")
@click.option("-n", "--n_samples", type=int, default=1, help="number of samples")
@click.option("-fc", "--force_constants", type=str, help="file with force constants")
@click.option("--mc_rattle", is_flag=True, help="use `mc_rattle` from hiphive")
@click.option("--quantum", is_flag=True, help="use quantum distribution function")
@click.option("--deterministic", is_flag=True, help="create a deterministic sample")
@click.option("-seed", "--random_seed", type=int, help="seed the random numbers")
@click.option("--format", default="aims")
def tool_create_samples(
    filename,
    temperature,
    n_samples,
    force_constants,
    mc_rattle,
    quantum,
    deterministic,
    random_seed,
    format,
):
    """create samples from geometry in FILENAME"""

    click.echo("hilde CLI: create_samples")
    create_samples(
        filename,
        temperature,
        n_samples,
        force_constants,
        mc_rattle,
        quantum,
        deterministic,
        random_seed,
        format,
    )


@tools.command("suggest_k_grid")
@click.argument("filename")
@click.option("-d", "--density", default=3.5)
@click.option("--uneven", is_flag=True)
@click.option("--format", default="aims")
def tool_suggest_k_grid(filename, density, uneven, format):
    """suggest a k_grid for geometry in FILENAME based on density"""

    click.echo("hilde CLI: suggest_k_grid")
    suggest_k_grid(filename, density, uneven, format)


@tools.command("remap_phonopy_force_constants")
@click.argument("filename")
@click.option("-uc", "--uc_filename", default="geometry.in.primitive")
@click.option("-sc", "--sc_filename", default="geometry.in.supercell")
def tool_remap_phonopy_force_constants(filename, uc_filename, sc_filename):
    """remap phonopy force constants in FILENAME to [3N, 3N] shape"""

    remap_phonopy_force_constants(
        uc_filename=uc_filename, sc_filename=sc_filename, fc_filename=filename
    )


@tools.command("nomad_upload")
@click.argument("folders", nargs=-1)
@click.option("--token", help="nomad token, otherwise read from .hilderc")
@click.option("--dry", help="only show the commands")
def tool_nomad_upload(folders, token, dry):
    """upload the calculations in FOLDERS to NOMAD"""

    nomad_upload(folders, token, dry)
