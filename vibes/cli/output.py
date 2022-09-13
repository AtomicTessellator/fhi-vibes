"""`vibes output` part of the CLI"""
from pathlib import Path

import click
import numpy as np

from vibes import defaults
from vibes import dimensions as dims
from vibes import keys
from vibes.filenames import filenames

from .misc import ClickAliasedGroup as AliasedGroup
from .misc import complete_files, default_context_settings


_default_context_settings = {"show_default": True}


@click.command(cls=AliasedGroup)
def output():
    """produce output of vibes workfow"""


@output.command(aliases=["md"], context_settings=default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-ha", "--harmonic", is_flag=True, help="comp. ha. properties incl. flux")
@click.option("-fc", "--fc_file", type=Path, help="add force constants from file")
@click.option("-o", "--outfile", default="auto", show_default=True)
@click.option("--force", is_flag=True, help="enfore parsing of output file")
@click.option("--shorten", default=0.0, help="shorten trajectory by percentage. Discard the first steps with positive value, the last with negative.")
def trajectory(file, harmonic, fc_file, outfile, force, shorten):
    """write trajectory data in FILE to xarray.Dataset"""
    if "auto" in outfile.lower():
        outfile = Path(file).stem
        outfile += ".nc"
    outfile = Path(outfile)

    file_size = Path(file).stat().st_size
    if not force and outfile.exists():
        import xarray as xr

        click.echo(f"Check if {file} has been parsed already")
        file_size_is = xr.open_dataset(outfile).attrs.get(keys.st_size)

        if file_size == file_size_is:
            click.echo(f".. input file (size) has not changed, skip.")
            click.echo(".. (use --force to parse anyway)")
            return
        else:
            click.echo(".. file size has changed, parse the file.")

    click.echo(f"Extract Trajectory dataset from {file}")
    from vibes.trajectory import reader
    from vibes.trajectory.dataset import get_trajectory_dataset
    from vibes import io

    traj = reader(file=file, fc_file=fc_file)
    if shorten > 0:
        click.echo(f".. shorten trajectory by {shorten*100} % from the first steps:")
        tmax_ds = float(traj.times[-1])
        click.echo(f"... max. time in trajectory: {tmax_ds} fs")
        n_max = len(traj)
        n_shorten = int(np.floor(n_max * shorten))
        click.echo(f"... discard first {n_shorten} steps")
        traj = traj.discard(first=n_shorten,last=0)
        traj.times = traj.times - traj.times[0]
        click.echo(f"... new trajectory length: {traj.times[-1]:.2f} fs")
        fc = io.parse_force_constants(fc_file, two_dim=False)
        traj.set_force_constants(fc)
    elif shorten < 0:
        shorten = -shorten
        click.echo(f".. shorten trajectory by {shorten*100} % from the last steps:")
        tmax_ds = float(traj.times[-1])
        click.echo(f"... max. time in trajectory: {tmax_ds} fs")
        n_max = len(traj)
        n_shorten = int(np.floor(n_max * shorten))
        click.echo(f"... discard last {n_shorten} steps")
        traj = traj.discard(first=0,last=n_shorten)
        click.echo(f"... new trajectory length: {traj.times[-1]:.2f} fs")
        fc = io.parse_force_constants(fc_file, two_dim=False)
        traj.set_force_constants(fc)

    if traj.stresses_potential is not None:
        traj.compute_heat_flux()
    if harmonic and traj.force_constants is not None:
        traj.compute_heat_flux_harmonic()

    DS = get_trajectory_dataset(traj, metadata=True)
    # attach file size
    DS.attrs.update({keys.st_size: file_size})
    # write to disk
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
    from vibes.phonopy import _defaults
    from vibes.phonopy.postprocess import extract_results, plot_results, postprocess

    phonon = postprocess(
        trajectory_file=file,
        born_charges_file=born,
        enforce_sum_rules=sum_rules,
    )
    if not q_mesh:
        if phonon.mesh_numbers is None:
            q_mesh = _defaults.kwargs.q_mesh.copy()
            click.echo(f"q_mesh not given, use default {q_mesh}")
        else:
            q_mesh = list(phonon.mesh_numbers)
            click.echo(f"q_mesh not given, use values stored in {file}: {q_mesh}")

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
    from vibes.phono3py.postprocess import extract_results, postprocess

    if not q_mesh:
        q_mesh = kwargs.q_mesh.copy()
        click.echo(f"q_mesh not given, use default {q_mesh}")

    phonon = postprocess(trajectory=file)

    output_directory = Path(file).parent / "output"

    extract_results(phonon, output_dir=output_directory)


@output.command(aliases=["gk"], context_settings=_default_context_settings)
@click.argument("file", default="trajectory.nc")
@click.option("-o", "--outfile", default="greenkubo.nc", type=Path)
@click.option("-w", "--window_factor", default=defaults.window_factor)
@click.option("--filter_prominence", default=defaults.filter_prominence)
@click.option("--interpolate", is_flag=True, help="interpolate to dense grid")
@click.option("--total", is_flag=True, help="compute total flux")
@click.option("-fc", "--fc_file", type=Path, help="use force constants from file")
@click.option("-u", "--update", is_flag=True, help="only parse if input data changed")
@click.option("--shorten", default=0.0, help="shorten trajectory by percentage. Discard the first steps with positive value, the last with negative.")
def greenkubo(
    file,
    outfile,
    window_factor,
    filter_prominence,
    interpolate,
    total,
    fc_file,
    update,
    shorten,
):
    """perform greenkubo analysis for dataset in FILE"""
    import xarray as xr

    import vibes.green_kubo as gk
    from vibes.io import parse_force_constants

    if total:
        outfile = outfile.parent / f"{outfile.stem}.total.nc"

    click.echo(f"Run aiGK output workflows for {file}")

    with xr.open_dataset(file) as ds_raw:

        ds = ds_raw
        if shorten > 0:
            dim = dims.time
            click.echo(f".. shorten trajectory by {shorten*100} % from the first steps:")
            tmax_ds = float(ds_raw[dim].isel({dim: -1}))
            click.echo(f"... max. time in trajectory: {tmax_ds} fs")
            n_max = len(ds_raw[dim])
            n_shorten = int(np.floor(n_max * shorten))
            click.echo(f"... discard first {n_shorten} steps")
            # ds = ds_raw.shift({dim: n_shorten}).dropna(dim=dim)  #ZS: will kill every step with NA
            ds = ds_raw.isel(time=slice(None,n_max-n_shorten))  #ZS: another shorten method
            ds = ds.assign_coords({dim: ds[dim] - ds[dim][0]})
            tm = float(ds.time.max() / 1000)
            click.echo(f"... new trajectory length: {tm*1000} fs")
        elif shorten < 0:
            shorten = -shorten
            dim = dims.time
            click.echo(f".. shorten trajectory by {shorten*100} % from the last steps:")
            tmax_ds = float(ds_raw[dim].isel({dim: -1}))
            click.echo(f"... max. time in trajectory: {tmax_ds} fs")
            n_max = len(ds_raw[dim])
            n_shorten = int(np.floor(n_max * shorten))
            click.echo(f"... discard last {n_shorten} steps")
            ds = ds_raw.isel(time=slice(n_shorten,None))
            ds = ds.assign_coords({dim: ds[dim] - ds[dim][0]})
            tm = float(ds.time.max() / 1000)
            click.echo(f"... new trajectory length: {tm*1000} fs")

        # check if postprocess is necessary
        if Path(outfile).exists() and update:
            file_size_old = xr.open_dataset(outfile).attrs.get(keys.st_size)
            file_size_new = ds.attrs.get(keys.st_size)

            if file_size_new == file_size_old:
                click.echo(f".. input file (size) has not changed, skip.")
                click.echo(".. (use --force to parse anyway)")
                return
            else:
                click.echo(f".. file size has changed, parse the file.")

        if fc_file is not None and keys.fc in ds:
            click.echo(f".. update force constants from {fc_file}")
            fcs = ds[keys.fc]
            fcs.data = parse_force_constants(fc_file)
            fcs.attrs = {"filename": str(Path(fc_file).absolute())}
            ds[keys.fc] = fcs

        ds_gk = gk.get_gk_dataset(
            ds,
            interpolate=interpolate,
            window_factor=window_factor,
            filter_prominence=filter_prominence,
            total=total,
        )

    click.echo(f".. write to {outfile}")

    ds_gk.to_netcdf(outfile)
