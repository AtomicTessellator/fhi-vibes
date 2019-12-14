"""`vibes output` part of the CLI"""
from pathlib import Path

# import numpy as np

import click

# from vibes.trajectory import reader
# from vibes.phonopy._defaults import defaults
# from vibes.tdep.wrapper import convert_phonopy_to_tdep
# from vibes.io import parse_force_constants

from .misc import AliasedGroup, complete_filenames


@click.command(cls=AliasedGroup)
def output():
    """produce output of vibes workfow"""


@output.command("md")
@click.argument("trajectory", default="trajectory.son", type=complete_filenames)
@click.option("-hf", "--heat_flux", is_flag=True, help="write heat flux dataset")
@click.option("-d", "--discard", type=int, help="discard this many steps")
@click.option("--minimal", is_flag=True, help="only write necessary minimum")
@click.option("-fc", "--force_constants", help="use FC to compute ha. forces")
@click.option("-rfc", "--remapped_force_constants", help="use remapped FC")
@click.option("-o", "--outfile", default="auto", show_default=True)
@click.option("-avg", "--average_reference", is_flag=True)
def md_output(
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
    import numpy as np
    from vibes.trajectory import reader
    from vibes.io import parse_force_constants

    click.echo(f"Extract Trajectory dataset from {trajectory}")
    traj = reader(file=trajectory, with_stresses=heat_flux, fc_file=force_constants)

    if discard:
        traj = traj.discard(discard)

    # harmonic forces?
    if force_constants:
        traj.set_forces_harmonic(average_reference=average_reference)
    elif remapped_force_constants:
        fc = np.loadtxt(remapped_force_constants)
        traj.set_force_constants_remapped(fc)
        traj.set_forces_harmonic(average_reference=average_reference)

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
    from vibes.tdep.wrapper import convert_phonopy_to_tdep

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

    # if tdep:
    #     convert_phonopy_to_tdep(phonon, workdir=str(output_directory) + "_tdep")
