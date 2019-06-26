"""`hilde output` part of the CLI"""
from pathlib import Path

import click

from hilde.phonopy._defaults import defaults

from .misc import AliasedGroup


@click.command(cls=AliasedGroup)
def output():
    """produce output of hilde workfow"""
    pass


@output.command("phonopy")
@click.argument("trajectory")
# necessary?
@click.option("--q_mesh", nargs=3, default=None)
@click.option("-od", "--output_directory")
@click.option("-bs", "--bandstructure", is_flag=True)
@click.option("-dos", "--density_of_states", is_flag=True)
@click.option("-tp", "--thermal_properties", is_flag=True)
@click.option("--animate", is_flag=True, help="print animation files for special kpts")
@click.option("--animate_q", nargs=3, multiple=True, type=float, help="animation at q")
@click.option("--full", is_flag=True)
@click.option("--tdep", is_flag=True, hidden=True)
@click.pass_obj
def phonopy_output(
    obj,
    trajectory,
    q_mesh,
    output_directory,
    bandstructure,
    density_of_states,
    thermal_properties,
    animate,
    animate_q,
    full,
    tdep,
):
    """perform phonopy postprocess for TRAJECTORY

    Parameters
    ----------
    obj: CliTracker
        The click context passed as an object
    trajectory: str
        Path to the phonopy trajectory file
    q_mesh: list of ints
        Size of the q-point interpolation mesh
    output_directory: str
        The output directory (default: output/)
    bandstructure: bool
        If True write and plot the bandstructure
    density_of_states: bool
        If True write and plot the total and projected density of states
    thermal_properties:
        If True write the thermal properties
    animate: bool
        If True write phonon mode animation files for all special q-points
    animate_q: list of tuples (float, float, float)
        A list of q-points to write an animation file for
    full: bool
        If True perform all output
    tdep: bool
        If True setup tdep results
    """
    from hilde.phonopy.postprocess import postprocess, extract_results

    if not q_mesh:
        q_mesh = defaults.q_mesh.copy()
        click.echo(f"q_mesh not given, use default {q_mesh}")

    phonon = postprocess(trajectory=trajectory)

    if not output_directory:
        output_directory = Path(trajectory).parent / "output"

    kwargs = {
        "write_thermal_properties": thermal_properties or full,
        "write_bandstructure": bandstructure or full,
        "write_dos": density_of_states or full,
        "write_pdos": density_of_states or full,
        "plot_bandstructure": bandstructure or full,
        "plot_dos": density_of_states or full,
        "plot_pdos": density_of_states or full,
        "q_mesh": q_mesh,
        "output_dir": output_directory,
        "tdep": tdep,
        "animate": animate or full,
        "animate_q": animate_q,
    }

    extract_results(phonon, **kwargs)
