"""`vibes output` part of the CLI"""
from pathlib import Path

import click

from vibes.filenames import filenames

from .misc import ClickAliasedGroup as AliasedGroup
from .misc import complete_files, default_context_settings


@click.command(cls=AliasedGroup)
def output():
    """produce output of vibes workfow"""


@output.command(aliases=["md"], context_settings=default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-hf", "--heat_flux", is_flag=True, help="write heat flux dataset")
@click.option("-d", "--discard", type=int, help="discard this many steps")
@click.option("--minimal", is_flag=True, help="omit redundant information")
@click.option("-fc", "--fc_file", type=Path, help="add force constants from file")
@click.option("-o", "--outfile", default="auto", show_default=True)
def trajectory(file, heat_flux, discard, minimal, fc_file, outfile):
    """write trajectory data in FILE to xarray.Dataset"""
    from vibes.trajectory import reader
    from vibes.trajectory.dataset import get_trajectory_dataset

    click.echo(f"Extract Trajectory dataset from {trajectory}")
    traj = reader(file=file, fc_file=fc_file)

    if discard:
        traj = traj.discard(discard)

    # harmonic forces?
    if fc_file:
        traj.set_forces_harmonic()

    if heat_flux:
        traj.compute_heat_fluxes_from_stresses()

    if "auto" in outfile.lower():
        outfile = Path(file).stem
        outfile += ".nc"

    DS = get_trajectory_dataset(traj, metadata=True)
    DS.to_netcdf(outfile)
    click.echo(f"Trajectory dataset written to {outfile}")


@output.command(context_settings=default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-bs", "--bandstructure", is_flag=True, help="plot bandstructure")
@click.option("--dos", is_flag=True, help="plot DOS")
@click.option("--full", is_flag=True, help="include thermal properties and animation")
@click.option("--q_mesh", nargs=3, default=None, help="use this q-mesh")
@click.option("--debye", is_flag=True, help="compute Debye temperature")
@click.option("-pdos", "--projected_dos", is_flag=True, help="plot projected DOS")
@click.option("--born", type=complete_files, help="include file with BORN charges")
@click.option("--sum_rules", is_flag=True, help="enfore sum rules with hiphive")
@click.option("-v", "--verbose", is_flag=True, help="print frequencies at gamma point")
@click.pass_obj
def phonopy(
    obj,
    file,
    bandstructure,
    dos,
    full,
    q_mesh,
    debye,
    projected_dos,
    born,
    sum_rules,
    verbose,
):
    """perform phonopy postprocess for trajectory in FILE"""
    from vibes.phonopy import _defaults as defaults
    from vibes.phonopy.postprocess import postprocess, extract_results, plot_results

    if not q_mesh:
        q_mesh = defaults.kwargs.q_mesh.copy()
        click.echo(f"q_mesh not given, use default {q_mesh}")

    phonon = postprocess(
        trajectory_file=file, born_charges_file=born, enforce_sum_rules=sum_rules,
    )

    folder = "output"
    if sum_rules:
        folder += "_sum_rules"
    output_directory = Path(file).parent / folder

    kwargs = {
        "minimal_output": True,
        "thermal_properties": full,
        "bandstructure": bandstructure or full,
        "dos": dos or full,
        "debye": debye,
        "pdos": projected_dos,
        "q_mesh": q_mesh,
        "output_dir": output_directory,
        "animate": full,
        "verbose": verbose,
    }

    extract_results(phonon, **kwargs)

    kwargs = {
        "thermal_properties": full,
        "bandstructure": bandstructure or full,
        "dos": dos or full,
        "pdos": projected_dos,
        "output_dir": output_directory,
    }
    plot_results(phonon, **kwargs)


@output.command()
@click.argument("file", default="trajectory.son", type=complete_files)
# necessary?
@click.option("--q_mesh", nargs=3, default=None)
@click.pass_obj
def phono3py(obj, file, q_mesh):
    """perform phono3py postprocess for trajectory in FILE"""
    from vibes.phono3py._defaults import kwargs
    from vibes.phono3py.postprocess import postprocess, extract_results

    if not q_mesh:
        q_mesh = kwargs.q_mesh.copy()
        click.echo(f"q_mesh not given, use default {q_mesh}")

    phonon = postprocess(trajectory=file)

    output_directory = Path(file).parent / "output"

    extract_results(phonon, output_dir=output_directory)


@output.command(aliases=["gk"])
@click.argument("file", default="trajectory_hf.nc")
@click.option("-avg", "--average", default=100, help="average window")
@click.option("--full", is_flag=True)
@click.option("--aux", is_flag=True)
@click.option("-o", "--outfile", default="greenkubo.nc", show_default=True, type=Path)
@click.option("-d", "--discard", default=0)
def greenkubo(file, average, full, aux, outfile, discard):
    """perform greenkubo analysis for dataset in FILE"""
    import xarray as xr
    import vibes.green_kubo.heat_flux as hf

    ds = xr.load_dataset(file)

    ds_kappa = hf.get_kappa_cumulative_dataset(ds, full=full, aux=aux, discard=discard)

    if full:
        outfile = outfile.parent / f"{outfile.stem}_full.nc"

    click.echo(f".. write to {outfile}")
    ds_kappa.to_netcdf(outfile)
