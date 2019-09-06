"""`hilde output` part of the CLI"""
from pathlib import Path

import click

from hilde.trajectory import reader
from hilde.phonopy._defaults import defaults
from hilde.io import parse_force_constants

from .misc import AliasedGroup, complete_filenames


@click.command(cls=AliasedGroup)
def output():
    """produce output of hilde workfow"""


@output.command("md")
@click.argument("trajectory", default="trajectory.son", type=complete_filenames)
@click.option("-hf", "--heat_flux", is_flag=True, help="write heat flux dataset")
@click.option("-d", "--discard", type=int, help="discard this many steps")
@click.option("--minimal", is_flag=True, help="only write necessary minimum")
@click.option("-fc", "--force_constants", help="use FC to compute ha. forces")
@click.option("-o", "--outfile", default="auto", show_default=True)
def md_output(trajectory, heat_flux, discard, minimal, force_constants, outfile):
    """write data in trajectory as xarray.Dataset"""

    click.echo(f"Extract Trajectory dataset from {trajectory}")
    traj = reader(file=trajectory, with_stresses=heat_flux)

    if discard:
        traj = traj.discard(discard)

    # harmonic forces?
    if force_constants:
        fc = parse_force_constants(
            force_constants, primitive=traj.primitive, supercell=traj.supercell
        )
        traj.set_forces_harmonic(force_constants=fc)

    if "auto" in outfile.lower():
        outfile = Path(trajectory).stem + ".nc"

    DS = traj.dataset
    DS.to_netcdf(outfile)
    click.echo(f"Trajectory dataset written to {outfile}")

    if heat_flux:
        outfile = "heat_flux.nc"
        DS = traj.get_heat_flux_data(only_flux=minimal)
        DS.to_netcdf(outfile)
        click.echo(f"Heat flux dataset written to {outfile}")


@output.command("phonopy")
@click.argument("trajectory", default="trajectory.son", type=complete_filenames)
# necessary?
@click.option("--q_mesh", nargs=3, default=None)
@click.option("-od", "--output_directory")
@click.option("-bs", "--bandstructure", is_flag=True)
@click.option("-dos", "--density_of_states", is_flag=True)
@click.option("-pdos", "--projected_density_of_states", is_flag=True)
@click.option("-tp", "--thermal_properties", is_flag=True)
@click.option("--animate", is_flag=True, help="print animation files for special kpts")
@click.option("--animate_q", nargs=3, multiple=True, type=float, help="animation at q")
@click.option("--born", type=complete_filenames)
@click.option("--full", is_flag=True)
@click.option("--tdep", is_flag=True, hidden=True)
@click.option("--remap_fc", is_flag=True, hidden=True)
@click.option("-v", "--verbose", is_flag=True, help="print frequencies at gamma point")
@click.pass_obj
def phonopy_output(
    obj,
    trajectory,
    q_mesh,
    output_directory,
    bandstructure,
    density_of_states,
    projected_density_of_states,
    thermal_properties,
    animate,
    animate_q,
    born,
    full,
    tdep,
    remap_fc,
    verbose,
):
    """perform phonopy postprocess for TRAJECTORY"""
    from hilde.phonopy.postprocess import postprocess, extract_results

    if not q_mesh:
        q_mesh = defaults.q_mesh.copy()
        click.echo(f"q_mesh not given, use default {q_mesh}")

    phonon = postprocess(
        trajectory=trajectory,
        born_charges_file=born,
        calculate_full_force_constants=remap_fc,
    )

    if not output_directory:
        output_directory = Path(trajectory).parent / "output"

    kwargs = {
        "write_thermal_properties": thermal_properties or full,
        "write_bandstructure": bandstructure or full,
        "write_dos": density_of_states or full,
        "write_pdos": projected_density_of_states,
        "plot_bandstructure": bandstructure or full,
        "plot_thermal_properties": thermal_properties or full,
        "plot_dos": density_of_states or full,
        "plot_pdos": projected_density_of_states,
        "q_mesh": q_mesh,
        "output_dir": output_directory,
        "tdep": tdep,
        "animate": animate or full,
        "animate_q": animate_q,
        "remap_fc": remap_fc,
        "verbose": verbose,
    }

    extract_results(phonon, **kwargs)
