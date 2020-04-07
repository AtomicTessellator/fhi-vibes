"""`vibes info` backend"""
from pathlib import Path

from vibes.filenames import filenames

from .misc import ClickAliasedGroup, click, complete_files


@click.command(cls=ClickAliasedGroup)
def info():
    """inform about content of a file"""


@info.command("settings")
@click.argument("file", type=complete_files)
@click.pass_obj
def settings_info(obj, file):
    """inform about content of a settings file"""
    from vibes.settings import Settings

    click.echo(f"List content of {file} including system-wide configuration")

    settings = Settings(file)
    settings.print()


@info.command("geometry")
@click.argument("file", default=filenames.atoms, type=complete_files)
@click.option("--format", default="aims", show_default=True)
@click.option("-t", "--symprec", default=1e-5, show_default=True)
@click.option("-v", "--verbose", is_flag=True)
@click.pass_obj
def geometry_info(obj, file, format, symprec, verbose):
    """inform about a structure in a geometry input file"""
    from ase.io import read
    from vibes.structure.io import inform

    atoms = read(file, format=format)

    verbosity = 1
    if verbose:
        verbosity = 2

    inform(atoms, symprec=symprec, verbosity=verbosity)


@info.command("md")
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-p", "--plot", is_flag=True, help="plot a summary")
@click.option("-w", "--write", is_flag=True, help="write Dataset to nc file")
@click.option("--avg", default=100, help="window size for running avg")
@click.option("-v", "--verbose", is_flag=True, help="be verbose")
def md_info(file, plot, write, avg, verbose):
    """inform about content of a settings.in file"""
    import xarray as xr
    from .scripts.md_sum import md_sum
    from vibes.trajectory import analysis as al, reader

    file = Path(file)

    if file.suffix in (".son", ".yaml", ".bz", ".gz"):
        trajectory = reader(file)
        DS = trajectory.dataset
        if write:
            trajectory.write(file.parent / f"{file.stem}.nc")
    elif file.suffix in (".nc"):
        DS = xr.load_dataset(file)
    elif file.suffix in (".log"):
        md_sum(file, plot, avg, verbose)
    else:
        raise click.FileError(f"File format of {file} not known.")

    click.echo(f"Dataset summary for {file}:")
    al.summary(DS, plot=plot, avg=avg)


@info.command("phonopy")
@click.argument("file", default="phonopy.in", type=complete_files)
@click.option("--write_supercell", is_flag=True, help="write the supercell to file")
@click.option("--format", default="aims", show_default=True)
def phonopy_info(file, write_supercell, format):
    """inform about a phonopy calculation before it is started"""
    from .scripts.vibes_phonopy import preprocess

    preprocess(
        file=None,
        settings_file=file,
        dimension=None,
        format=format,
        write_supercell=write_supercell,
    )


@info.command("trajectory")
@click.argument("file", default=filenames.trajectory, type=complete_files)
def trajectory_info(file):
    """inform about content of trajectory file"""
    from vibes import son
    from vibes.settings import Settings

    metadata, _ = son.load(file)

    click.echo(f"Summary of metadata in {file}:\n")
    click.echo("Keys:")
    click.echo(f"  {list(metadata.keys())}\n")
    if "settings" in metadata:
        settings = Settings.from_dict(metadata["settings"])
        click.echo("Settings:")
        settings.print()


@info.command("netcdf")
@click.argument("file", type=complete_files)
def show_netcdf_file(file):
    """show contents of netCDF FILE"""
    import xarray as xr

    DS = xr.open_dataset(file)

    print(DS)


@info.command("hdf5")
@click.argument("file", type=complete_files)
@click.option("-v", "--verbose", is_flag=True)
def show_hdf5_file(file, verbose):
    """show contents of HDF5 FILE"""
    import pandas as pd

    click.echo(f"Summarize file {file}")
    with pd.HDFStore(file) as store:
        click.echo("Keys:")
        for k in store:
            click.echo(f"  {k}")

        if verbose:
            click.echo()
            for k in store:
                click.echo(f"Describe {k}")
                click.echo(store[k].describe())


@info.command(aliases=["gk"])
@click.argument("dataset", default="greenkubo.nc")
@click.option("-p", "--plot", is_flag=True, help="plot summary")
@click.option("--no_hann", is_flag=True)
@click.option("--logx", is_flag=True)
@click.option("--xlim", type=float, help="xlim range in ps")
@click.option("-avg", "--average", default=100, help="average window")
def greenkubo(dataset, plot, no_hann, logx, xlim, average):
    import xarray as xr
    from vibes import keys
    from vibes.green_kubo.analysis import summary, plot_summary

    DS = xr.load_dataset(dataset)

    (df_time, df_freq) = summary(DS)

    if plot:
        fig = plot_summary(
            df_time,
            df_freq,
            t_avalanche=DS.attrs[keys.time_avalanche],
            logx=logx,
            xlim=xlim,
            avg=average,
        )

        file = Path(dataset).stem + "_summary.pdf"
        fig.savefig(file, bbox_inches="tight")
        click.echo(f".. summary plotted to {file}")


@info.command("vdos")
@click.argument("file", default=filenames.trajectory_dataset, type=complete_files)
@click.option("-o", "--output_file", default="vdos.csv")
@click.option("-p", "--plot", is_flag=True, help="plot the DOS")
@click.option("--peak", type=float, help="height for peak detection", show_default=1)
@click.option("-mf", "--max_frequency", default=30.0, help="max. freq. in THz")
def velocity_autocorrelation(file, output_file, plot, peak, max_frequency):
    """write velocity autocorrelation function to output file"""
    import xarray as xr
    from vibes.green_kubo.velocities import get_vdos, simple_plot

    click.echo(f"Read {file} and extract velocities")
    velocities = xr.open_dataset(file).velocities

    vdos = get_vdos(velocities=velocities, hann=False, verbose=True)

    # sum atoms and coordinates
    df = vdos.real.sum(axis=(1, 2)).to_series()

    if plot:
        simple_plot(df, height=peak, max_frequency=max_frequency)

    click.echo(f".. write VDOS to {output_file}")
    df.to_csv(output_file, index_label="omega", header=True)


@info.command("relaxation")
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-v", "--verbose", is_flag=True)
@click.pass_obj
def relaxation_info(obj, file, verbose):
    """inform about a structure in a geometry input file"""
    from ase.constraints import full_3x3_to_voigt_6_stress
    from vibes.relaxation.context import MyExpCellFilter as ExpCellFilter
    from vibes.trajectory import reader
    from vibes.spglib.wrapper import get_spacegroup
    from vibes.relaxation._defaults import keys, kwargs, relaxation_options, name

    traj, metadata = reader(file, get_metadata=True, verbose=False)

    relaxation_kwargs = metadata[name].get(relaxation_options, {})

    try:
        fmax = relaxation_kwargs[keys.fmax]
    except KeyError:
        fmax = kwargs[keys.fmax]
        click.echo(f"** fmax not found in {file}, use default value {fmax}")

    try:
        fix_symmetry = relaxation_kwargs[keys.fix_symmetry]
    except KeyError:
        fix_symmetry = kwargs[keys.fix_symmetry]
        msg = f"** `fix_symmetry` not found in {file}, use default value {fix_symmetry}"
        click.echo(msg)

    scalar_pressure = relaxation_kwargs.get(
        keys.scalar_pressure, kwargs[keys.scalar_pressure]
    )

    atoms_ref = traj[0]
    na = len(atoms_ref)

    energy_ref = atoms_ref.get_potential_energy()

    click.echo(f"Relaxation info for {file}:")
    if verbose:
        import json

        click.echo("Metadata for relaxation:")
        click.echo(json.dumps(metadata[name], indent=2))

    if fix_symmetry:
        from ase.spacegroup.symmetrize import FixSymmetry

        click.echo("fix_symmetry:     True")

    if scalar_pressure:
        click.echo(f"scalar_pressure: {scalar_pressure*1000: .3e} meV/A**3")

    click.echo(f"fmax:            {fmax*1000: .3e} meV/AA")
    click.echo(
        "# Step |   Free energy   |   F-F(1)   | max. force |  max. stress |"
        + "  Volume  |  Spacegroup  |"
        + "\n"
        + "#      |       [eV]      |    [meV]   |  [meV/AA]  |  [meV/AA^3]  |"
        + "  [AA^3]  |              |"
        + "\n"
    )

    for ii, atoms in enumerate(traj[1:]):

        energy = atoms.get_potential_energy()
        de = 1000 * (energy - energy_ref)

        opt_atoms = ExpCellFilter(atoms, scalar_pressure=scalar_pressure)

        forces = opt_atoms.get_forces()
        stress = full_3x3_to_voigt_6_stress(forces[na:])
        forces = forces[:na]  # drop the stress

        # optionally: symmetrize forces and stress
        if fix_symmetry:
            constr = FixSymmetry(atoms, symprec=kwargs[keys.symprec])
            constr.adjust_forces(atoms, forces)
            constr.adjust_stress(atoms, stress)

        res_forces = (forces ** 2).sum(axis=1).max() ** 0.5 * 1000
        res_stress = abs(stress).max() * 1000

        vol_str = f"{atoms.get_volume():10.3f}"
        # sg_str = f"{get_spacegroup(atoms):5d}"
        sg_str = get_spacegroup(atoms)

        msg = "{:5d}   {:16.8f}  {:12.6f} {:12.4f} {:14.4f} {}   {}".format(
            ii + 1, energy, de, res_forces, res_stress, vol_str, sg_str,
        )
        click.echo(msg)

    if max(res_forces, res_stress) < fmax * 1000:
        click.echo("--> converged.")
