"""`vibes output` part of the CLI"""
from pathlib import Path

import click

from vibes.filenames import filenames

from .misc import ClickAliasedGroup as AliasedGroup
from .misc import complete_files


@click.command(cls=AliasedGroup)
def output():
    """produce output of vibes workfow"""


@output.command(aliases=["md"])
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-hf", "--heat_flux", is_flag=True, help="write heat flux dataset")
@click.option("-d", "--discard", type=int, help="discard this many steps")
@click.option("--minimal", is_flag=True, help="only write necessary minimum")
@click.option("-fc", "--force_constants", help="use FC to compute ha. forces")
@click.option("-o", "--outfile", default="auto", show_default=True)
@click.option("-avg", "--average_reference", is_flag=True)
def trajectory(
    file, heat_flux, discard, minimal, force_constants, outfile, average_reference,
):
    """write trajectory data in FILE to xarray.Dataset"""
    from vibes.trajectory import reader
    from vibes.trajectory.dataset import get_trajectory_dataset

    click.echo(f"Extract Trajectory dataset from {trajectory}")
    traj = reader(file=file, fc_file=force_constants)

    if discard:
        traj = traj.discard(discard)

    # harmonic forces?
    if force_constants:
        traj.set_forces_harmonic(average_reference=average_reference)

    if heat_flux:
        traj.compute_heat_fluxes_from_stresses()

    if "auto" in outfile.lower():
        outfile = Path(file).stem
        outfile += ".nc"

    DS = get_trajectory_dataset(traj, metadata=True)
    DS.to_netcdf(outfile)
    click.echo(f"Trajectory dataset written to {outfile}")


@output.command()
@click.argument("file", default=filenames.trajectory, type=complete_files)
# necessary?
@click.option("--q_mesh", nargs=3, default=None)
@click.option("-bs", "--bandstructure", is_flag=True)
@click.option("-dos", "--density_of_states", is_flag=True)
@click.option("-debye", "--debye_temperature", is_flag=True)
@click.option("-pdos", "--projected_density_of_states", is_flag=True)
@click.option("-tp", "--thermal_properties", is_flag=True)
@click.option("-path", "--bz_path", type=str)
@click.option("--animate", is_flag=True, help="print animation files for special kpts")
@click.option("--animate_q", nargs=3, multiple=True, type=float, help="animation at q")
@click.option("--born", type=complete_files)
@click.option("--full", is_flag=True)
@click.option("--remap_fc", is_flag=True)
@click.option("--sum_rules", is_flag=True)
@click.option("-v", "--verbose", is_flag=True, help="print frequencies at gamma point")
@click.pass_obj
def phonopy(
    obj,
    file,
    q_mesh,
    bandstructure,
    density_of_states,
    debye_temperature,
    projected_density_of_states,
    thermal_properties,
    bz_path,
    animate,
    animate_q,
    born,
    full,
    remap_fc,
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
        trajectory_file=file,
        born_charges_file=born,
        calculate_full_force_constants=remap_fc,
        enforce_sum_rules=sum_rules,
    )

    folder = "output"
    if sum_rules:
        folder += "_sum_rules"
    output_directory = Path(file).parent / folder

    kwargs = {
        "minimal_output": True,
        "thermal_properties": thermal_properties or full,
        "bandstructure": bandstructure or full,
        "dos": density_of_states or full,
        "debye": debye_temperature or full,
        "pdos": projected_density_of_states,
        "bz_path": bz_path,
        "q_mesh": q_mesh,
        "output_dir": output_directory,
        "animate": animate or full,
        "animate_q": animate_q,
        "remap_fc": remap_fc,
        "verbose": verbose,
    }

    extract_results(phonon, **kwargs)

    kwargs = {
        "thermal_properties": thermal_properties or full,
        "bandstructure": bandstructure or full,
        "dos": density_of_states or full,
        "pdos": projected_density_of_states,
        "bz_path": bz_path,
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
