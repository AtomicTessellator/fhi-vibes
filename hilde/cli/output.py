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
@click.option("--full", is_flag=True)
@click.option("--tdep", is_flag=True, hidden=True)
@click.option("--animate", is_flag=True, help="print animation files for special kpts")
@click.option("--animate_q", nargs=3, multiple=True, type=float, help="animation at q")
@click.pass_obj
def phonopy_output(
    obj,
    trajectory,
    q_mesh,
    output_directory,
    bandstructure,
    density_of_states,
    thermal_properties,
    full,
    tdep,
    animate,
    animate_q,
):
    """perform phonopy postprocess for TRAJECTORY"""
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
        "animate": animate,
        "animate_q": animate_q,
    }

    extract_results(phonon, **kwargs)
