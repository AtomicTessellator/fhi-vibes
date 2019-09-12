"""hilde CLI utils"""

from pathlib import Path
from hilde.scripts.refine_geometry import refine_geometry
from hilde.scripts.make_supercell import make_supercell
from hilde.scripts.create_samples import create_samples
from hilde.scripts.suggest_k_grid import suggest_k_grid
from hilde.scripts.remap_phonopy_forceconstants import remap_force_constants
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
@click.option("-d", "--dimension", type=int, nargs=9)
@click.option("-dd", "--diagonal_dimension", type=int, nargs=3)
@click.option("-n", "--n_target", type=int)
@click.option("--deviation", default=0.1, show_default=True)
@click.option("--dry", is_flag=True)
@click.option("--format", default="aims")
@click.option("--scaled", is_flag=True)
def tool_make_supercell(
    filename, dimension, diagonal_dimension, n_target, deviation, dry, format, scaled
):
    """create a supercell of desired shape or size"""
    if diagonal_dimension:
        dimension = diagonal_dimension
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


@utils.command("remap_force_constants")
@click.argument("filename", default="FORCE_CONSTANTS", type=complete_filenames)
@click.option("-uc", "--primitive", default="geometry.in.primitive", show_default=True)
@click.option("-sc", "--supercell", default="geometry.in.supercell", show_default=True)
@click.option("--python", is_flag=True)
def tool_remap_phonopy_force_constants(filename, primitive, supercell, python):
    """remap phonopy force constants in FILENAME to [3N, 3N] shape"""

    remap_force_constants(
        fc_filename=filename,
        uc_filename=primitive,
        sc_filename=supercell,
        fortran=not python,
    )


@utils.command("nomad_upload")
@click.argument("folders", nargs=-1, type=complete_filenames)
@click.option("--token", help="nomad token, otherwise read from .hilderc")
@click.option("--legacy", is_flag=True, default=True, help="use old Nomad (default)")
@click.option("--dry", is_flag=True, help="only show the commands")
def tool_nomad_upload(folders, token, legacy, dry):
    """upload the calculations in FOLDERS to NOMAD"""

    nomad_upload(folders, token, legacy, dry)


@utils.command(cls=AliasedGroup)
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
@click.option("-uc", help="Add a (primitive) unit cell", type=complete_filenames)
@click.option("-sc", help="Add the respective supercell", type=complete_filenames)
@click.option("--format", default="aims")
def trajectory_update(filename, uc, sc, format):
    """add unit cell from UC and supercell from SC to trajectory in FILENAME"""
    update_trajectory(filename, uc, sc, format)


@trajectory.command("pick_sample")
@click.argument("filename", default="trajectory.son", type=complete_filenames)
@click.option("-n", "--number", default=0)
@click.option("-cart", "--cartesian", is_flag=True, help="write cart. coords")
def pick_sample(filename, number, cartesian):
    """pick a sample from trajectory and write to geometry input file"""
    from hilde.trajectory import reader

    click.echo(f"Read trajectory from {filename} and extract sample {number}:")

    traj = reader(filename)

    if number < 0:
        number = len(traj) + number

    atoms = traj[number]

    outfile = f"geometry.in.{number}"
    atoms.write(outfile, format="aims", velocities=True, scaled=not cartesian)

    click.echo(f".. sample written to {outfile}")


@utils.command(cls=AliasedGroup)
def harmonicity():
    """utils for quantifying harmonicity"""
    ...


@harmonicity.command("r2")
@click.argument("filenames", type=complete_filenames, nargs=-1)
@click.option("-o", "--outfile", help="Save dataframe to csv")
@click.option("--quiet", is_flag=True)
@click.option("--pick", type=int, help="pick one sample")
def compute_r2(filenames, outfile, quiet, pick):
    """Compute r2 and some statistics"""

    from pathlib import Path
    import pandas as pd
    import xarray as xr
    from hilde.anharmonicity_score import get_r, get_r2_per_sample, get_r2_MAE

    click.echo(f"Compute harmonicity score for {len(filenames)} materials:")

    dct = {}
    for filename in filenames:
        click.echo(f" parse {filename}")

        DS = xr.open_dataset(filename)

        f = DS.forces
        fh = DS.forces_harmonic

        if pick is not None:
            f = DS.forces[pick]
            fh = DS.forces_harmonic[pick]
        else:
            f = DS.forces
            fh = DS.forces_harmonic

        r = get_r(f, fh)
        r2s = get_r2_per_sample(f, fh)
        r2a = get_r2_MAE(f, fh)

        name = DS.attrs["System Name"]
        dct[name] = {
            "r": r,
            "r2 (RMS)": r2s.mean(),
            "r2 (RMA)": 1 - (1 - r2a) ** 2,
            "fa (MSE)": (1 - r2s.mean()) ** 0.5,
            "fa (MAE)": 1 - r2a,
        }

    df = pd.DataFrame(dct).T
    if not quiet:
        click.echo("\nDataFrame:")
        click.echo(df)
        click.echo("\nDataFrame.describe():")
        click.echo(df.describe())

    if outfile is not None:
        inp = f"{outfile} exists, overwrite? (Y/n) "
        if Path(outfile).exists() and input(inp) == "Y":
            return
        df.to_csv(outfile, index_label="material", float_format="%15.12e")
        click.echo(f"\n.. Dataframe written to {outfile}")


@utils.command(cls=AliasedGroup)
def hdf5():
    """utils for working with hdf5 files and pandas"""
    ...


@hdf5.command("describe")
@click.argument("file", type=complete_filenames)
@click.option("-v", "--verbose", is_flag=True)
def show_hdf5_file(file, verbose):
    """show contents of HDF5 FILE"""
    import pandas as pd

    click.echo(f"Summarize file {file}")
    with pd.HDFStore(file) as store:
        click.echo("Keys:")
        for k in store:
            click.echo(f"  {k}")

        if verbose:
            click.echo()
            for k in store:
                click.echo(f"Describe {k}")
                click.echo(store[k].describe())


@hdf5.command("join")
@click.argument("files", nargs=-1, type=complete_filenames)
@click.option("-o", "--outfile", default="data.h5", show_default=True)
def join_hdf5(files, outfile):
    """join csv containing pandas.DataFrames to one hdf5 store"""
    import pandas as pd

    click.echo(f"Join files {files}")

    with pd.HDFStore(outfile) as store:
        for file in files:
            click.echo(f'.. append {file} to {outfile}')
            df = pd.read_csv(file)
            key = Path(file).stem.replace(".", "_")
            store[key] = df
