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


def dict2atoms(atoms_dict):
    from ase.atoms import Atoms
    from ase.symbols import symbols2numbers
    import json

    atoms_dict = json.loads(atoms_dict)
    atoms_dict["numbers"] = symbols2numbers(atoms_dict["symbols"])
    new_atoms_dict = {}
    for index in ["numbers", "positions", "cell"]:
        new_atoms_dict[index] = np.array(atoms_dict[index])
    new_atoms_dict["pbc"] = atoms_dict["pbc"]
    return Atoms.fromdict(new_atoms_dict)


@click.command(cls=AliasedGroup)
def output():
    """produce output of vibes workfow"""


@output.command(aliases=["md"], context_settings=default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-ha", "--harmonic", is_flag=True, help="comp. ha. properties incl. flux")
@click.option("-fc", "--fc_file", type=Path, help="add force constants from file")
@click.option("-o", "--outfile", default="auto", show_default=True)
@click.option("--force", is_flag=True, help="enfore parsing of output file")
@click.option("--shorten", default=0.0, help="shorten trajectory by percentage")
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

    traj = reader(file=file)

    if shorten != 0:
        click.echo(f".. shorten trajectory by {shorten*100} %")
        tmax_ds = float(traj.times[-1])
        n_max = len(traj)
        n_shorten = int(np.floor(n_max * abs(shorten)))
        if shorten > 0:
            traj = traj.discard(first=n_shorten, last=0)
            traj.times = traj.times - traj.times[0]
        elif shorten < 0:
            traj = traj.discard(first=0, last=n_shorten)
        click.echo(f"... max. time in trajectory: {tmax_ds} fs")
        click.echo(f"... new trajectory length: {traj.times[-1]:.2f} fs")
    if fc_file:
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
@click.option("-c", "--cross_offdiag", is_flag=True, help="cross offdiag cutoff time")
@click.option("-fc", "--fc_file", type=Path, help="use force constants from file")
@click.option("--fc_supercell", type=Path, help="interpolate with new FC & supercell")
@click.option("--born", type=complete_files, help="include file with BORN charges")
@click.option("-u", "--update", is_flag=True, help="only parse if input data changed")
@click.option("--shorten", default=0.0, help="shorten trajectory by percentage")
def greenkubo(
    file,
    outfile,
    window_factor,
    filter_prominence,
    interpolate,
    total,
    cross_offdiag,
    fc_file,
    fc_supercell,
    born,
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

    with xr.open_dataset(file) as ds:
        if shorten != 0:
            dim = dims.time
            click.echo(f".. shorten trajectory by {shorten*100} %")
            tmax_ds = float(ds[dim].isel({dim: -1}))
            click.echo(f"... max. time in trajectory: {tmax_ds} fs")
            n_max = len(ds[dim])
            n_shorten = int(np.floor(n_max * abs(shorten)))
            if shorten > 0:
                click.echo(f"... discard first {n_shorten} steps")
                ds = ds.isel(time=slice(n_shorten, None))
            elif shorten < 0:
                click.echo(f"... discard last {n_shorten} steps")
                ds = ds.isel(time=slice(None, n_max - n_shorten))
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
            if fc_supercell is not None:
                from ase.io import read
                from vibes.helpers.converters import atoms2json
                import vibes.dimensions as dims

                click.echo(f".. use another force constant for interpolation")
                interpolation_supercell = atoms2json(
                    read(fc_supercell, format="aims"), reduce=False
                )
                interpolation_fc = parse_force_constants(fc_file)
                fc_dims = (dims.i, "J1", dims.a, dims.b)
                fcs = xr.DataArray(interpolation_fc, dims=fc_dims)
                fcs.attrs = {"filename": str(Path(fc_file).absolute())}
                ds[keys.interpolation_fc] = fcs
                ds.attrs[keys.interpolation_supercell] = interpolation_supercell
                click.echo(f".. interpolation force constants from {fc_file}")
                click.echo(f".. interpolation supercell from {fc_supercell}")
            else:
                click.echo(f".. update force constants from {fc_file}")
                fcs = ds[keys.fc]
                fcs.data = parse_force_constants(fc_file)
                fcs.attrs = {"filename": str(Path(fc_file).absolute())}
                ds[keys.fc] = fcs

        if born is not None:
            click.echo(f".. update BEC from {born}")
            ds.attrs["born_charges_file"] = born
            click.echo(f".. update BEC from {ds.attrs['born_charges_file']}")

        ds_gk = gk.get_gk_dataset(
            ds,
            interpolate=interpolate,
            window_factor=window_factor,
            filter_prominence=filter_prominence,
            cross_offdiag=cross_offdiag,
            total=total,
        )

    click.echo(f".. write to {outfile}")

    ds_gk.to_netcdf(outfile)
