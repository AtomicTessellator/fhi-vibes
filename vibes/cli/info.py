"""`vibes info` backend"""
import warnings
from pathlib import Path

import numpy as np

from vibes import defaults, keys
from vibes.filenames import filenames

from .misc import ClickAliasedGroup, click, complete_files


# from click 7.1 on
_default_context_settings = {"show_default": True}


@click.command(cls=ClickAliasedGroup)
def info():
    """inform about content of a file"""


@info.command()
@click.argument("file", type=complete_files)
@click.pass_obj
def settings(obj, file):
    """write the settings in FILE *including* the configuration"""
    from vibes.settings import Settings

    click.echo(f"List content of {file} including system-wide configuration")

    settings = Settings(file)
    settings.print()


@info.command()
@click.argument("file", default=filenames.atoms, type=complete_files)
@click.option("--format", default="aims", show_default=True)
@click.option("-t", "--symprec", default=1e-5, show_default=True)
@click.option("-v", "--verbose", is_flag=True, help="increase verbosity")
@click.pass_obj
def geometry(obj, file, format, symprec, verbose):
    """inform about a structure in a geometry input file"""
    from ase.io import read

    from vibes.structure.io import inform

    atoms = read(file, format=format)

    verbosity = 1
    if verbose:
        verbosity = 2

    inform(atoms, symprec=symprec, verbosity=verbosity)


@info.command(context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory_dataset, type=complete_files)
@click.option("-p", "--plot", is_flag=True, help="plot a summary")
@click.option("--avg", default=100, help="window size for running avg")
@click.option("-v", "--verbose", is_flag=True, help="be verbose")
def md(file, plot, avg, verbose):
    """inform about MD simulation in FILE"""
    import xarray as xr

    from vibes.trajectory import analysis as al
    from vibes.trajectory import reader

    from .scripts.md_sum import md_sum

    file = Path(file)

    if file.suffix in (".son", ".yaml", ".bz", ".gz"):
        trajectory = reader(file)
        DS = trajectory.dataset
    elif file.suffix in (".nc"):
        DS = xr.load_dataset(file)
    elif file.suffix in (".log"):
        md_sum(file, plot, avg, verbose)
    else:
        raise click.FileError(f"File format of {file} not known.")

    click.echo(f"Dataset summary for {file}:")
    al.summary(DS, plot=plot, avg=avg)


@info.command(context_settings=_default_context_settings)
@click.argument("file", default="phonopy.in", type=complete_files)
@click.option("--write_supercell", is_flag=True, help="write the supercell to file")
def phonopy(file, write_supercell):
    """inform about a phonopy calculation based on the input FILE"""
    from .scripts.vibes_phonopy import preprocess

    preprocess(settings_file=file, write_supercell=write_supercell)


@info.command()
@click.argument("file", default=filenames.trajectory, type=complete_files)
def trajectory(file):
    """print metadata from trajectory in FILE"""
    from vibes import son
    from vibes.settings import Settings

    metadata, _ = son.open(file)

    click.echo(f"Summary of metadata in {file}:\n")
    click.echo("Keys:")
    click.echo(f"  {list(metadata.keys())}\n")
    if "settings" in metadata:
        settings = Settings.from_dict(metadata["settings"])
        click.echo("Settings:")
        settings.print()


@info.command()
@click.argument("file", type=complete_files)
@click.option("--key", type=str, default=None, help="print DS[key]")
def netcdf(file, key):
    """show contents of netCDF FILE"""
    import xarray as xr

    DS = xr.open_dataset(file)

    if key is not None:
        print(DS[key])
    else:
        print(DS)


@info.command(context_settings=_default_context_settings)
@click.argument("file", type=complete_files)
@click.option("--describe", is_flag=True, help="print description of data")
@click.option("--half", is_flag=True, help="print only the second half of data")
@click.option("--to_json", type=Path, help="Write to json file")
def csv(file, describe, half, to_json):
    """show contents of csv FILE"""
    import json

    import pandas as pd

    df = pd.read_csv(file)

    pd.options.display.max_rows = len(df)
    pd.options.display.max_columns = 0  # automatic detection

    if half:
        df = df.iloc[len(df) // 2 :]

    if describe:
        df = df.describe()

    click.echo(df)

    if to_json is not None:
        click.echo(f".. write to {to_json}")
        with open(to_json, "w") as f:
            json.dump(df.to_dict(), f, indent=1)


@info.command(aliases=["gk"], context_settings=_default_context_settings)
@click.argument("files", nargs=-1)
@click.option("-p", "--plot", is_flag=True, help="plot info")
@click.option("--lifetimes", is_flag=True, help="describe lifetimes")
# @click.option("--png", is_flag=True, help="plot png instead of pdf")
def greenkubo(files, plot, lifetimes):
    """visualize heat flux and thermal conductivity"""
    import xarray as xr

    if len(files) < 1:
        files = ["greenkubo.nc"]

    datasets = [xr.load_dataset(f) for f in files]

    # concatentate trajectories
    ds_gk = xr.concat(datasets, dim=keys.trajectory)
    # ds_gk = ds_gk.sel({keys.time: datasets[0][keys.time]})
    attrs = {ii: d.attrs for (ii, d) in enumerate(datasets)}
    ds_gk.attrs = attrs

    ks = ds_gk[keys.kappa]
    ks_flat = ks.stack(ab=("a", "b"))[:, ::4].data
    k = ks_flat.mean()
    k_err = (ks_flat.var() / (ks_flat.size)) ** 0.5
    click.echo(f"Kappa:    {k:.3f} +/- {k_err:.3f}")
    click.echo(f"Kappa^ab: {ks.mean(axis=0)}")

    # harmonic/interpolation
    if keys.kappa_ha in ds_gk:
        if lifetimes:
            click.echo("Lifetimes:")
            click.echo(ds_gk[keys.mode_lifetime].to_series().describe())

        click.echo("Harmonic values:")
        ks = ds_gk[keys.kappa_ha]
        ks_flat = ks.stack(ab=("a", "b"))[:, ::4].data
        k_ha = ks_flat.mean()
        k_ha_err = (ks_flat.var() / (ks_flat.size)) ** 0.5
        click.echo(f"Kappa:    {k_ha:.3f} +/- {k_ha_err:.3f}")
        click.echo(f"Kappa^ab: {ks.mean(axis=0)}")

        if len(datasets) > 1:
            click.echo(".. perform ensemble average before proceding")

        # ignore warnings regarding empty slices
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ds_gk = ds_gk.mean(dim=keys.trajectory)

        # interpolation
        intercept = float(ds_gk[keys.interpolation_fit_intercept])
        correction1 = float(ds_gk[keys.interpolation_correction])
        correction2 = intercept - k_ha

        click.echo("Summarize interpolation:")
        click.echo(f"Initial harmonic kappa value:    {k_ha:.3f} W/mK")
        click.echo(f"Fit intercept:                   {intercept:.3f} W/mK")
        click.echo(f"Interpolated harm. kappa:        {k_ha+correction1:.3f} W/mK")
        click.echo(f"Correction1 (k_init + m / nq):   {correction1:.3f} W/mK")
        click.echo(f"Correction2 (k_inf  - k_init):   {correction2:.3f} W/mK")

    if plot:
        from .scripts.plot_gk_summary import main as plot_summary

        outfile = Path(files[0]).stem + "_summary.pdf"
        plot_summary(ds_gk, outfile=outfile)

        if keys.interpolation_correction in ds_gk:
            from .scripts.plot_gk_interpolation import main as plot_interpolation

            _outfile = Path(outfile).stem + "_interpolation.pdf"
            plot_interpolation(ds_gk, outfile=_outfile)


@info.command(context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory_dataset, type=complete_files)
@click.option("-o", "--outfile", default="vdos.csv")
@click.option("-p", "--plot", is_flag=True, help="plot the DOS")
@click.option("--filter_prominence", default=defaults.filter_prominence)
@click.option("-mf", "--max_frequency", type=float, help="max. freq. in THz")
@click.option("--npad", default=10000, help="number of zeros for padding")
def vdos(file, outfile, plot, filter_prominence, max_frequency, npad):
    """compute and write velocity autocorrelation function to output file"""
    import xarray as xr

    from vibes.green_kubo.velocities import get_vdos, simple_plot

    click.echo(f"Read {file} and extract velocities")
    velocities = xr.open_dataset(file).velocities

    vdos = get_vdos(velocities=velocities, hann=False, npad=npad, verbose=True)

    # sum atoms and coordinates
    df = vdos.real.sum(axis=(1, 2)).to_series()

    if plot:
        peaks = simple_plot(
            df, prominence=filter_prominence, max_frequency=max_frequency
        )
        file = Path(outfile).stem + "_peaks.dat"
        click.echo(f".. save peak positions to {file}")
        np.savetxt(file, peaks)

    click.echo(f".. write VDOS to {outfile}")
    df.to_csv(outfile, index_label="omega", header=True)


@info.command(context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-v", "--verbose", is_flag=True, help="show more information")
@click.pass_obj
def relaxation(obj, file, verbose):
    """summarize geometry optimization in FILE"""
    from ase.constraints import full_3x3_to_voigt_6_stress

    from vibes.relaxation._defaults import keys, kwargs, name, relaxation_options
    from vibes.relaxation.context import MyExpCellFilter as ExpCellFilter
    from vibes.spglib import get_spacegroup
    from vibes.trajectory import reader

    traj = reader(file, verbose=False)

    relaxation_kwargs = traj.metadata[name].get(relaxation_options, {})

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
        click.echo(json.dumps(traj.metadata[name], indent=2))

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


@info.command(context_settings=_default_context_settings)
@click.argument("files", nargs=-1, type=complete_files)
@click.option("--per_sample", is_flag=True, help="analyze per sample")
@click.option("--per_mode", is_flag=True, help="analyze per mode (no --per_sample)")
@click.option("--describe", is_flag=True)
@click.option("--dry", is_flag=True, help="don't write output files")
def anharmonicity(files, per_sample, per_mode, describe, dry):
    """Compute sigmaA for trajectory dataset in FILE"""
    import pandas as pd
    import xarray as xr

    from vibes import keys
    from vibes.anharmonicity_score import get_dataframe, get_sigma_per_mode

    if len(files) < 1:
        files = [filenames.trajectory_dataset]

    dfs = []
    for file in files:
        click.echo(f"Compute anharmonicity measure for {file}:")

        click.echo(f" parse {file}")

        DS = xr.open_dataset(file)

        outfile = "sigmaA"

        name = DS.attrs[keys.system_name]
        if per_mode:
            df = get_sigma_per_mode(DS)
            outfile += "_mode"
            index_label = keys.omega
        else:
            df = get_dataframe(DS, per_sample=per_sample)
            index_label = "material"
        dfs.append(df)

    df = pd.concat(dfs)

    click.echo("\nDataFrame:")
    click.echo(df)
    if describe:
        click.echo("\nDataFrame.describe():")
        click.echo(df.describe())

    if outfile is not None and not dry:
        df.to_csv(f"{outfile}.csv", index_label=index_label, float_format="%15.12e")
        click.echo(f"\n.. Dataframe for {name} written to {outfile}")
