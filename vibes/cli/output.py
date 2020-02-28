"""`vibes output` part of the CLI"""
from pathlib import Path

import click

from .misc import ClickAliasedGroup as AliasedGroup
from .misc import complete_filenames


@click.command(cls=AliasedGroup)
def output():
    """produce output of vibes workfow"""


@output.command(aliases=["md"])
@click.argument("trajectory", default="trajectory.son", type=complete_filenames)
@click.option("-hf", "--heat_flux", is_flag=True, help="write heat flux dataset")
@click.option("-d", "--discard", type=int, help="discard this many steps")
@click.option("--minimal", is_flag=True, help="only write necessary minimum")
@click.option("-fc", "--force_constants", help="use FC to compute ha. forces")
@click.option("-rfc", "--remapped_force_constants", help="use remapped FC")
@click.option("-o", "--outfile", default="auto", show_default=True)
@click.option("-avg", "--average_reference", is_flag=True)
def molecular_dynamics(
    trajectory,
    heat_flux,
    discard,
    minimal,
    force_constants,
    remapped_force_constants,
    outfile,
    average_reference,
):
    """write data in trajectory as xarray.Dataset"""
    from vibes.trajectory import reader
    from vibes.io import parse_force_constants
    from vibes.trajectory.dataset import get_trajectory_dataset

    click.echo(f"Extract Trajectory dataset from {trajectory}")
    traj = reader(file=trajectory, fc_file=force_constants)

    if discard:
        traj = traj.discard(discard)

    # harmonic forces?
    if force_constants:
        traj.set_forces_harmonic(average_reference=average_reference)
    elif remapped_force_constants:
        fc = parse_force_constants(remapped_force_constants)
        traj.set_force_constants_remapped(fc)
        traj.set_forces_harmonic(average_reference=average_reference)

    if heat_flux:
        traj.compute_heat_fluxes_from_stresses()

    if "auto" in outfile.lower():
        file = Path(trajectory).stem
        if heat_flux:
            file += "_hf"
        outfile = file + ".nc"

    DS = get_trajectory_dataset(traj, metadata=True)
    DS.to_netcdf(outfile)
    click.echo(f"Trajectory dataset written to {outfile}")


@output.command("phonopy")
@click.argument("trajectory", default="trajectory.son", type=complete_filenames)
# necessary?
@click.option("--q_mesh", nargs=3, default=None)
@click.option("-od", "--output_directory")
@click.option("-bs", "--bandstructure", is_flag=True)
@click.option("-dos", "--density_of_states", is_flag=True)
@click.option("-debye", "--debye_temperature", is_flag=True)
@click.option("-pdos", "--projected_density_of_states", is_flag=True)
@click.option("-tp", "--thermal_properties", is_flag=True)
@click.option("-path", "--bz_path", type=str)
@click.option("--animate", is_flag=True, help="print animation files for special kpts")
@click.option("--animate_q", nargs=3, multiple=True, type=float, help="animation at q")
@click.option("--born", type=complete_filenames)
@click.option("--full", is_flag=True)
@click.option("--remap_fc", is_flag=True)
@click.option("--sum_rules", is_flag=True)
@click.option("-v", "--verbose", is_flag=True, help="print frequencies at gamma point")
@click.pass_obj
def phonopy_output(
    obj,
    trajectory,
    q_mesh,
    output_directory,
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
    """perform phonopy postprocess for TRAJECTORY"""
    from vibes.phonopy._defaults import defaults
    from vibes.phonopy.postprocess import postprocess, extract_results

    if not q_mesh:
        q_mesh = defaults.q_mesh.copy()
        click.echo(f"q_mesh not given, use default {q_mesh}")

    phonon = postprocess(
        trajectory=trajectory,
        born_charges_file=born,
        calculate_full_force_constants=remap_fc,
        enforce_sum_rules=sum_rules,
    )

    if not output_directory:
        folder = "output"
        if sum_rules:
            folder += "_sum_rules"
        output_directory = Path(trajectory).parent / folder

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


@output.command(aliases=["gk"])
@click.argument("dataset", default="trajectory_hf.nc")
@click.option("-avg", "--average", default=100, help="average window")
@click.option("--full", is_flag=True)
@click.option("--aux", is_flag=True)
@click.option("-o", "--outfile", default="greenkubo.nc", show_default=True, type=Path)
@click.option("-d", "--discard", default=0)
def greenkubo(dataset, average, full, aux, outfile, discard):
    """perform greenkubo analysis"""
    import xarray as xr
    import vibes.green_kubo.heat_flux as hf

    ds = xr.load_dataset(dataset)

    ds_kappa = hf.get_kappa_cumulative_dataset(ds, full=full, aux=aux, discard=discard)

    if full:
        outfile = outfile.parent / f"{outfile.stem}_full.nc"

    click.echo(f".. write to {outfile}")
    ds_kappa.to_netcdf(outfile)
