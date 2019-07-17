"""hilde CLI utils"""

from hilde.scripts.refine_geometry import refine_geometry
from hilde.scripts.make_supercell import make_supercell
from hilde.scripts.create_samples import create_samples
from hilde.scripts.suggest_k_grid import suggest_k_grid
from hilde.scripts.remap_phonopy_forceconstants import remap_phonopy_force_constants
from hilde.scripts.nomad_upload import nomad_upload
from hilde.scripts.update_md_trajectory import update_trajectory
from hilde.scripts.get_relaxation_info import get_relaxation_info

from .misc import click, AliasedGroup, complete_filenames


@click.command(cls=AliasedGroup)
def utils():
    """tools and utils"""


@utils.command(cls=AliasedGroup)
def geometry():
    """utils for working with structures"""
    ...


@geometry.command("refine")
@click.argument("filename", type=complete_filenames)
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


@utils.command("make_supercell")
@click.argument("filename", default="geometry.in", type=complete_filenames)
@click.option("-d", "--dimension", type=float, nargs=9)
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


@utils.command("get_relaxation_info")
@click.argument("filenames", nargs=-1, type=complete_filenames)
def relaxation_info(filenames):
    """analyze aims relaxation"""
    get_relaxation_info(filenames)


@utils.command("create_samples")
@click.argument("filename", type=complete_filenames)
@click.option("-T", "--temperature", type=float, help="Temperature in Kelvin")
@click.option("-n", "--n_samples", type=int, default=1, help="number of samples")
@click.option("-fc", "--force_constants", type=complete_filenames)
@click.option("--mc_rattle", is_flag=True, help="`hiphive.mc_rattle`", hidden=True)
@click.option("--quantum", is_flag=True, help="use quantum distribution function")
@click.option("--deterministic", is_flag=True, help="create a deterministic sample")
@click.option("--sobol", is_flag=True, help="use Sobol numbers to create samples")
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
    sobol,
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
        sobol,
        random_seed,
        format,
    )


@utils.command("suggest_k_grid")
@click.argument("filename", type=complete_filenames)
@click.option("-d", "--density", default=3.5)
@click.option("--uneven", is_flag=True)
@click.option("--format", default="aims")
def tool_suggest_k_grid(filename, density, uneven, format):
    """suggest a k_grid for geometry in FILENAME based on density"""

    click.echo("hilde CLI: suggest_k_grid")
    suggest_k_grid(filename, density, uneven, format)


@utils.command("remap_phonopy_force_constants")
@click.argument("filename", default="FORCE_CONSTANTS", type=complete_filenames)
@click.option("-pc", "--primitive", default="geometry.in.primitive", show_default=True)
@click.option("-sc", "--supercell", default="geometry.in.supercell", show_default=True)
def tool_remap_phonopy_force_constants(filename, primitive, supercell):
    """remap phonopy force constants in FILENAME to [3N, 3N] shape"""

    remap_phonopy_force_constants(
        uc_filename=primitive, sc_filename=supercell, fc_filename=filename
    )


@utils.command("nomad_upload")
@click.argument("folders", nargs=-1, type=complete_filenames)
@click.option("--token", help="nomad token, otherwise read from .hilderc")
@click.option("--dry", is_flag=True, help="only show the commands")
def tool_nomad_upload(folders, token, dry):
    """upload the calculations in FOLDERS to NOMAD"""

    nomad_upload(folders, token, dry)


@utils.command(cls=AliasedGroup, hidden=True)
def trajectory():
    """trajectory utils"""


@trajectory.command("2tdep")
@click.argument("filename", default="trajectory.son", type=complete_filenames)
@click.option("-s", "--skip", default=1, help="skip this many steps from trajectory")
@click.option("--folder", default="tdep", help="folder to store input")
def t2tdep(filename, skip, folder):
    """extract trajectory in FILENAME and store tdep input files to FOLDER"""
    from hilde.trajectory import reader

    traj = reader(filename)
    traj.to_tdep(folder=folder, skip=skip)


@trajectory.command("2xyz")
@click.argument("filename", default="trajectory.son", type=complete_filenames)
@click.option("--file", default="trajectory.xyz")
def t2xyz(filename, file):
    """extract trajectory in FILENAME and store as xyz file"""
    from hilde.trajectory import reader

    traj = reader(filename)
    traj.to_xyz(file=file)


@trajectory.command("update")
@click.argument("filename", default="trajectory.son", type=complete_filenames)
@click.option("-uc", help="Add a (primitive) unit cell")
@click.option("-sc", help="Add the respective supercell")
@click.option("--format", default="aims")
def trajectory_update(filename, uc, sc, format):
    """add unit cell from UC and supercell from SC to trajectory in FILENAME"""
    update_trajectory(filename, uc, sc, format)


@trajectory.command("extract_velocities")
@click.argument("filename", default="trajectory.son", type=complete_filenames)
@click.option("-o", "--output_filename", default="velocities.nc")
def extract_velocities(filename, output_filename):
    """extract velocities from FILENAME and write as xarray.DataArray to netCDF file"""
    from hilde.green_kubo.velocities import get_velocities

    velocities = get_velocities(trajectory=filename)

    velocities.to_netcdf(output_filename)
    click.echo(f".. velocities written to {output_filename}")
