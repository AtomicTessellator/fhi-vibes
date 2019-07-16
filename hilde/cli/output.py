"""`hilde output` part of the CLI"""
from pathlib import Path

import click

from hilde.phonopy._defaults import defaults

from .misc import AliasedGroup, complete_filenames


@click.command(cls=AliasedGroup)
def output():
    """produce output of hilde workfow"""


@output.command("phonopy")
@click.argument("trajectory", type=complete_filenames)
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
    verbose,
):
    """perform phonopy postprocess for TRAJECTORY"""
    from hilde.phonopy.postprocess import postprocess, extract_results

    if not q_mesh:
        q_mesh = defaults.q_mesh.copy()
        click.echo(f"q_mesh not given, use default {q_mesh}")

    phonon = postprocess(trajectory=trajectory, born_charges_file=born)

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
        "verbose": verbose,
    }

    extract_results(phonon, **kwargs)
