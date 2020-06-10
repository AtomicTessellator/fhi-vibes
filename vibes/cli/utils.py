"""vibes CLI utils"""

from vibes import keys
from vibes.filenames import filenames

from .misc import AliasedGroup, ClickAliasedGroup, Path, click, complete_files

xrange = range

# from click 7.1 on
_default_context_settings = {"show_default": True}


@click.group(cls=ClickAliasedGroup)
def utils():
    """tools and utils"""


@utils.command(aliases=["hash"], context_settings=_default_context_settings)
@click.argument("file", type=complete_files)
@click.option("--dry", is_flag=True, help="Only print hash to stdout")
@click.option("-o", "--outfile", default="hash.toml")
def hash_file(file, dry, outfile):
    """create sha hash for FILE"""
    import time
    from vibes.helpers.utils import talk
    from vibes.helpers.hash import hashfunc

    timestr = time.strftime("%Y/%m/%d_%H:%M:%S")

    hash = hashfunc(open(file).read())
    talk(f'Hash for "{file}":\n  {hash}')

    if not dry:
        with open(outfile, "a") as f:
            f.write(f"\n# {timestr}")
            f.write(f'\n"{file}": "{hash}"')
        talk(f".. written to {outfile}")


@utils.command(cls=AliasedGroup)
def geometry():
    """utils for working with structures"""
    ...


@geometry.command(context_settings=_default_context_settings)  # aliases=['gd'])
@click.argument("files", nargs=2, type=complete_files)
@click.option("-sc_file", "--supercell_file", type=Path)
@click.option("--outfile", default=filenames.deformation, type=Path)
@click.option("--dry", is_flag=True)
@click.option("--format", default="aims")
def get_deformation(files, supercell_file, outfile, dry, format):
    """get matrix D that deformes lattice of geometry1 to lattice of geometry2, where
        A_2 = D A_1

    For Supercells:
        S = M A  =>  S_2 = M D A_1  = M D M^-1 S_1  = D_s S_1
    """
    import numpy as np
    from ase.io import read
    from vibes.helpers.geometry import get_deformation

    atoms1 = read(files[0], format=format)
    atoms2 = read(files[1], format=format)

    D = get_deformation(atoms1.cell, atoms2.cell).round(14)

    # supercell?
    if supercell_file is not None:
        click.echo(f"Reference supercell from: {supercell_file}")
        sc = read(supercell_file, format=format)
        M = get_deformation(atoms1.cell, sc.cell).round(14)

        D = M @ D @ np.linalg.inv(M)

        outfile = f"{outfile.stem}_supercell{outfile.suffix}"

    click.echo("Deformation tensor:")
    click.echo(D)

    if not dry:
        import numpy as np

        click.echo(f".. write to {outfile}")
        np.savetxt(outfile, D)


@geometry.command(context_settings=_default_context_settings)  # aliases=['ad'])
@click.argument("file", type=complete_files)
@click.option("--deformation", type=complete_files, default=filenames.deformation)
@click.option("--outfile", type=Path)
@click.option("--dry", is_flag=True)
@click.option("--cartesian", is_flag=True)
@click.option("--format", default="aims")
def apply_deformation(file, deformation, outfile, dry, cartesian, format):
    """get matrix that deformes geometry1 to geometry2"""
    import numpy as np
    from ase.io import read

    D = np.loadtxt(deformation)
    atoms = read(file, format=format)

    click.echo(f"Apply deformation tensor from {deformation}:")
    click.echo(D)

    new_atoms = atoms.copy()
    new_atoms.set_cell(D @ atoms.cell[:], scale_atoms=True)

    velocities = False
    if atoms.get_velocities() is not None:
        velocities = True
        new_atoms.set_velocities(atoms.get_velocities())
        click.echo("** velocities are not scaled")

    info_str = f"Deformation tensor: {D.tolist()}"

    if not dry:
        outfile = outfile or file + ".deformed"
        click.echo(f".. write structure to {outfile}")
        kw = {
            "info_str": info_str,
            "scaled": not cartesian,
            "format": format,
            "velocities": velocities,
        }
        new_atoms.write(outfile, **kw)


@geometry.command(context_settings=_default_context_settings)
@click.argument("files", nargs=2, type=complete_files)
@click.option("--cartesian", is_flag=True)
@click.option("--format", default="aims")
def set_cell(files, cartesian, format):
    """set cell of geometry in file1 to geometry in file2"""
    from ase.io import read

    atoms1 = read(files[0], format=format)
    atoms2 = read(files[1], format=format)

    velocities = atoms2.get_velocities()

    click.echo(f"Replace lattice in {files[1]} by lattice from {files[0]}")
    atoms2.set_cell(atoms1.cell[:], scale_atoms=True)
    if velocities is not None:
        click.echo("** velocities are not scaled")

    click.echo(f".. write structure to {files[1]}")
    kw = {
        "scaled": not cartesian,
        "format": format,
        "velocities": velocities is not None,
    }
    atoms2.write(files[1], **kw)


@geometry.command("2frac", context_settings=_default_context_settings)
@click.argument("file", default=filenames.atoms, type=complete_files)
@click.option("-o", "--output_file")
@click.option("--format", default="aims")
def to_fractional(file, output_file, format):
    """rewrite geometry in fractional coordinates"""
    from ase.io import read

    if not output_file:
        output_file = file + ".fractional"

    atoms = read(file, format=format)
    atoms.write(output_file, format=format, scaled=True, geo_constrain=True)

    click.echo(f"Geometry written to {output_file}")


@geometry.command("2cart", context_settings=_default_context_settings)
@click.argument("file", default=filenames.atoms, type=complete_files)
@click.option("-o", "--output_file")
@click.option("--format", default="aims")
def to_cartesian(file, output_file, format):
    """rewrite geometry in cartesian coordinates"""
    from ase.io import read

    if not output_file:
        output_file = file + ".cartesian"

    atoms = read(file, format=format)
    atoms.write(output_file, format=format, scaled=False)

    click.echo(f"Geometry written to {output_file}")


@geometry.command("refine", context_settings=_default_context_settings)
@click.argument("file", type=complete_files)
@click.option("-prim", "--primitive", is_flag=True)
@click.option("-conv", "--conventional", is_flag=True)
@click.option("--center", is_flag=True)
@click.option("--origin", is_flag=True)
@click.option("-cart", "--cartesian", is_flag=True)
@click.option("--format", default="aims")
@click.option("-t", "--symprec", default=1e-5)
def geometry_refine(*args, **kwargs):
    """vibes.scripts.refine_geometry"""
    from .scripts.refine_geometry import refine_geometry

    refine_geometry(*args, **kwargs)


@geometry.command("wrap", context_settings=_default_context_settings)
@click.argument("file", default=filenames.atoms, type=complete_files)
@click.option("-o", "--output_file")
@click.option("--format", default="aims")
def wrap_atoms(file, output_file, format):
    """rewrite geometry in fractional coordinates"""
    from ase.io import read

    if not output_file:
        output_file = file + ".wrapped"

    atoms = read(file, format=format)
    atoms.positions -= [0.1, 0.1, 0.1]
    atoms.wrap(pretty_translation=True)
    atoms.positions += [0.1, 0.1, 0.1]
    atoms.wrap()
    atoms.write(output_file, format=format, scaled=True, geo_constrain=True, wrap=False)

    click.echo(f"Wrapped geometry written to {output_file}")


@utils.command("make_supercell", context_settings=_default_context_settings)
@click.argument("file", default=filenames.atoms, type=complete_files)
@click.option("-d", "--dimension", type=int, nargs=9)
@click.option("-dd", "--diagonal_dimension", type=int, nargs=3)
@click.option("-n", "--n_target", type=int)
@click.option("-o", "--output_file")
@click.option("--deviation", default=0.2)
@click.option("--dry", is_flag=True)
@click.option("--format", default="aims")
@click.option("-frac", "--fractional", is_flag=True)
@click.option("--wrap", is_flag=False)
def tool_make_supercell(
    file,
    dimension,
    diagonal_dimension,
    output_file,
    n_target,
    deviation,
    dry,
    format,
    fractional,
    wrap,
):
    """create a supercell of desired shape or size"""
    from .scripts.make_supercell import make_supercell

    if diagonal_dimension:
        dimension = diagonal_dimension

    make_supercell(
        file,
        dimension,
        n_target,
        deviation,
        dry,
        format,
        fractional,
        output_file=output_file,
        wrap=wrap,
    )


@utils.group()
def aims():
    """utils for working with FHI-aims (output)"""
    ...


@aims.command(context_settings=_default_context_settings)
@click.argument("files", nargs=-1, type=complete_files)
def get_relaxation_info(files):
    """analyze aims relaxation"""
    from .scripts.get_relaxation_info import get_relaxation_info

    get_relaxation_info(files)


@utils.command(aliases=["samples"], context_settings=_default_context_settings)
@click.argument("filename", type=complete_files)
@click.option("-T", "--temperature", type=float, help="Temperature in Kelvin")
@click.option("-n", "--n_samples", type=int, default=1, help="number of samples")
@click.option("-fc", "--force_constants_file", type=complete_files)
@click.option("--rattle", type=float, help="atoms.rattle(stdev=X) (ASE default: 0.001)")
@click.option("--quantum", is_flag=True, help="use quantum distribution function")
@click.option("--deterministic", is_flag=True, help="create a deterministic sample")
@click.option("--gauge_eigenvectors", is_flag=True)
@click.option("--zacharias", is_flag=True, help="Zacharias Sampling (deterministic)")
@click.option("--ignore_negative", is_flag=True)
@click.option("-seed", "--random_seed", type=int, help="seed the random numbers")
@click.option("--propagate", type=float, help="propagate this many fs")
@click.option("--format", default="aims")
def create_samples(filename, **kwargs):
    """create samples from geometry in FILENAME"""
    from .scripts.create_samples import create_samples

    click.echo("vibes CLI: create_samples")
    create_samples(atoms_file=filename, **kwargs)


@utils.command("suggest_k_grid", context_settings=_default_context_settings)
@click.argument("file", type=complete_files)
@click.option("-d", "--density", default=3.5)
@click.option("--uneven", is_flag=True)
@click.option("--format", default="aims")
def tool_suggest_k_grid(file, density, uneven, format):
    """suggest a k_grid for geometry in FILENAME based on density"""
    from .scripts.suggest_k_grid import suggest_k_grid

    click.echo("vibes CLI: suggest_k_grid")
    suggest_k_grid(file, density, uneven, format)


@utils.group(aliases=["fc"])
def force_constants():
    """utils for working with force constants"""
    ...


@force_constants.command(context_settings=_default_context_settings)
@click.argument("file", default="FORCE_CONSTANTS", type=complete_files)
@click.option("-uc", "--primitive", default=filenames.primitive)
@click.option("-sc", "--supercell", default=filenames.supercell)
@click.option("-nsc", "--new_supercell")
@click.option("-o", "--output_file")
@click.option("--symmetrize", is_flag=True)
@click.option("--python", is_flag=True)
@click.option("--format", default="aims")
def remap(
    file,
    primitive,
    supercell,
    new_supercell,
    output_file,
    symmetrize,
    python,
    eps=1e-13,
    tol=1e-5,
    format="aims",
):
    """remap phonopy force constants in FILENAME to [3N, 3N] shape"""
    # copy: from vibes.scripts.remap_phonopy_forceconstants import remap_force_constants
    import numpy as np
    from ase.io import read
    from vibes.io import parse_force_constants
    from vibes.phonopy.utils import remap_force_constants

    uc = read(primitive, format=format)
    sc = read(supercell, format=format)

    nsc = None
    if new_supercell:
        nsc = read(new_supercell, format=format)

    kwargs = {
        "primitive": uc,
        "supercell": sc,
        "fortran": not python,
        "eps": eps,
        "tol": tol,
    }

    fc = parse_force_constants(fc_file=file, two_dim=False, **kwargs)

    kwargs.update({"new_supercell": nsc, "two_dim": True, "symmetrize": symmetrize})

    fc = remap_force_constants(fc, **kwargs)

    if not output_file:
        output_file = f"{file}_remapped"

    msg = f"remapped force constants from {file}, shape [{fc.shape}]"
    np.savetxt(output_file, fc, header=msg)

    click.echo(f".. remapped force constants written to {output_file}")


@force_constants.command(context_settings=_default_context_settings)
@click.argument("file", default="FORCE_CONSTANTS_remapped", type=complete_files)
@click.option("-sc", "--supercell", default=filenames.supercell)
@click.option("-n", "--show_n_frequencies", default=6, type=int)
@click.option("-o", "--output_file", default="frequencies.dat")
@click.option("--symmetrize", is_flag=True)
@click.option("--format", default="aims")
def frequencies(file, supercell, show_n_frequencies, output_file, symmetrize, format):
    """compute the frequency spectrum"""
    import numpy as np
    from ase.io import read
    from vibes.harmonic_analysis.dynamical_matrix import get_frequencies

    atoms = read(supercell, format=format)
    fc = np.loadtxt(file)

    w2 = get_frequencies(fc, masses=atoms.get_masses())

    if show_n_frequencies:
        nn = show_n_frequencies
        print(f"The first {nn} frequencies:")
        for ii, freq in enumerate(w2[:nn]):
            print(f" {ii + 1:4d}: {freq}")

        print(f"Highest {nn} frequencies")
        for ii, freq in enumerate(w2[-nn:]):
            print(f" {len(w2) - ii:4d}: {freq }")

    if isinstance(output_file, str):
        np.savetxt(output_file, w2)
        click.echo(f".. frequencies written to {output_file}")


@utils.group()
def nomad():
    ...


@nomad.command("upload", context_settings=_default_context_settings)
@click.argument("files", nargs=-1, type=complete_files)
@click.option("--token", help="nomad token, otherwise read from .vibesrc")
@click.option("--name", help="nomad upload name")
@click.option("--legacy", is_flag=True, help="use old Nomad")
@click.option("--dry", is_flag=True, help="only show the commands")
def tool_nomad_upload(files, token, name, legacy, dry):
    """upload FILES to NOMAD"""
    from .scripts.nomad_upload import nomad_upload

    nomad_upload(files, token, legacy, dry, name=name)


@utils.group(cls=ClickAliasedGroup)
def trajectory():
    """trajectory utils"""


@trajectory.command(context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-n", "--n_steps", type=int, help="discard this many steps")
def discard(file, n_steps):
    """discard this many steps from trajectory"""
    from vibes.trajectory import reader

    traj = reader(file)

    click.echo(f"Discard {n_steps} steps")
    traj = traj[n_steps:-1]

    p = Path(file)
    outfile = f"{p.stem}_discarded{p.suffix}"

    traj.write(outfile)


@trajectory.command("2tdep", context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-s", "--skip", default=1, help="skip this many steps from trajectory")
@click.option("--folder", default="tdep", help="folder to store input")
def t2tdep(file, skip, folder):
    """extract trajectory in FILENAME and store tdep input files to FOLDER"""
    from vibes.trajectory import reader

    traj = reader(file)
    traj.to_tdep(folder=folder, skip=skip)


@trajectory.command("2xyz", context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-o", "--output_file", default="trajectory.xyz")
def t2xyz(file, output_file):
    """extract trajectory in FILENAME and store as xyz file"""
    from vibes.trajectory import reader

    traj = reader(file)
    traj.to_xyz(file=output_file)


@trajectory.command("2db", context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-o", "--output_file", default="trajectory.db")
def t2db(file, output_file):
    """extract trajectory in FILENAME and store as ase db"""
    from vibes.trajectory import reader

    traj = reader(file)
    traj.to_db(output_file)


@trajectory.command("2csv", context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-o", "--output_file", default="trajectory.csv")
def t2csv(file, output_file):
    """extract trajectory in FILENAME and store 1D data as csv dataframe"""
    from vibes.trajectory import reader

    traj = reader(file)
    df = traj.dataframe

    click.echo(f"Write trajectory data to {output_file}")
    df.to_csv(output_file)


@trajectory.command(context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-uc", help="Add a (primitive) unit cell", type=complete_files)
@click.option("-sc", help="Add the respective supercell", type=complete_files)
@click.option("-fc", help="Add the force constants", type=complete_files)
@click.option("-o", "--output_file")
@click.option("--format", default="aims")
def update(file, uc, sc, fc, output_file, format):
    """add unit cell from UC and supercell from SC to trajectory in FILENAME"""
    # copy: from vibes.scripts.update_md_trajectory import update_trajectory
    import shutil
    from ase.io import read
    from vibes.trajectory import reader

    traj = reader(file, fc_file=fc)

    if uc:
        atoms = read(uc, format=format)
        traj.primitive = atoms

    if sc:
        atoms = read(sc, format=format)
        traj.supercell = atoms

    if not output_file:
        suffix = Path(file).suffix
        new_trajectory = f"temp{suffix}"
        fname = f"{file}.bak"
        click.echo(f".. back up old trajectory to {fname}")
        shutil.copy(file, fname)

    else:
        new_trajectory = output_file

    traj.write(file=new_trajectory)

    if not output_file:
        click.echo(f".. move new trajectory to {file}")
        shutil.move(new_trajectory, file)


@trajectory.command(aliases=["pick"], context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-o", "--outfile", type=Path)
@click.option("-n", "--number", default=0)
@click.option("-r", "--range", type=int, nargs=3, help="start, stop, step")
@click.option("-cart", "--cartesian", is_flag=True, help="write cart. coords")
def pick_samples(file, outfile, number, range, cartesian):
    """pick a sample from trajectory and write to geometry input file"""
    from vibes.trajectory import reader

    click.echo(f"Read trajectory from {file}:")
    traj = reader(file)

    if number < 0:
        number = len(traj) + number

    if len(range) == 3:
        rge = xrange(*range)
    else:
        rge = [number]

    for number in rge:
        click.echo(f"Extract sample {number}:")
        file = outfile or f"{filenames.atoms}.{number}"
        atoms = traj[number]
        info_str = f"Sample no.: {number:7d}"
        kw = {"velocities": True, "scaled": not cartesian, "info_str": info_str}
        atoms.write(file, format="aims", **kw)
        click.echo(f".. sample written to {file}")


@trajectory.command("average", context_settings=_default_context_settings)
@click.argument("file", default=filenames.trajectory, type=complete_files)
@click.option("-r", "--range", type=int, nargs=3, help="start, stop, step")
@click.option("-cart", "--cartesian", is_flag=True, help="write cart. coords")
@click.option("-o", "--outfile", default="geometry.in.average")
def average_trajectory(file, range, cartesian, outfile):
    """average positions"""
    from vibes.trajectory import reader

    click.echo(f"Read trajectory from {file}:")
    traj = reader(file)

    if len(range) == 3:
        rge = xrange(*range)
        traj = traj[rge]

    atoms = traj.average_atoms
    atoms.wrap(pretty_translation=True)
    atoms.write(outfile, format="aims", scaled=not cartesian)
    click.echo(f".. geometry with averaged positions written to {outfile}")


@utils.group(aliases=["a"])
def anharmonicity():
    """utils for quantifying anharmonicity"""
    ...


@anharmonicity.command(context_settings=_default_context_settings)
@click.argument("file", type=complete_files)
@click.option("-o", "--outfile")
@click.option("--per_sample", is_flag=True)
@click.option("--describe", is_flag=True)
def sigma(file, outfile, per_sample, describe):
    """Compute sigmaA for trajectory dataset in FILE"""
    import xarray as xr
    from vibes import keys
    from vibes.anharmonicity_score import get_dataframe

    click.echo(f"Compute anharmonicity measure for {file}:")

    click.echo(f" parse {file}")

    DS = xr.open_dataset(file)

    name = DS.attrs[keys.system_name]
    df = get_dataframe(DS, per_sample=per_sample,)

    click.echo("\nDataFrame:")
    click.echo(df)
    if describe:
        click.echo("\nDataFrame.describe():")
        click.echo(df.describe())

    if outfile is not None:
        df.to_csv(outfile, index_label="material", float_format="%15.12e")
        click.echo(f"\n.. Dataframe for {name} written to {outfile}")

    # if store_hdf5:
    #     outfile += ".h5"
    #     with pd.HDFStore(outfile) as store:
    #         click.echo(f".. append dataframe for {name} to {outfile}")
    #         store[name] = df


@anharmonicity.command(context_settings=_default_context_settings)
@click.argument("file", type=complete_files)
@click.option("-o", "--outfile")
@click.option("--describe", is_flag=True)
def mode(file, outfile, describe):
    """Compute sigmaA per mode for trajectory dataset in FILE"""
    import xarray as xr
    from vibes import keys
    from vibes.anharmonicity_score import get_sigma_per_mode

    click.echo(f"Compute anharmonicity measure per mode for {file}")

    DS = xr.open_dataset(file)

    name = DS.attrs[keys.system_name]

    series = get_sigma_per_mode(DS)

    if describe:
        click.echo("\nSeries.describe():")
        click.echo(series.describe())

    outfile = outfile or f"sigmaA_mode_{name}.csv"
    series.to_csv(outfile, index_label=keys.omega)
    click.echo(f"\n.. mode resolved sigmaA for {name} written to {outfile}")


@utils.command("backup", context_settings=_default_context_settings)
@click.argument("folder", type=complete_files)
@click.option("--target", default=keys.default_backup_folder)
@click.option("--nozip", is_flag=True)
def perform_backup(folder, target, nozip):
    """backup FOLDER to TARGET"""
    from vibes.helpers.backup import backup_folder

    backup_folder(folder, target_folder=target, zip=not nozip)


@utils.group(aliases=["ph3"])
def phono3py():
    """utils for working with phono3py"""
    ...


@phono3py.command(
    "run_thermal_conductivity", context_settings=_default_context_settings
)
@click.argument("folder", default="output", type=complete_files)
@click.option("--q_mesh", nargs=3, help="q_mesh")
@click.option("--outfile", default="kappa_QMESH.log")
@click.pass_obj
def run_thermal_conductivity(obj, folder, q_mesh, outfile):
    """run a phono3py thermal conductivity calculation in FOLDER"""
    import sys
    from vibes.helpers.utils import talk
    from vibes.phono3py._defaults import kwargs
    from vibes.cli.scripts import run_thermal_conductivity as rtc

    if q_mesh is None:
        q_mesh = kwargs.q_mesh

    outfile = outfile.replace("QMESH", ".".join((str(q) for q in q_mesh)))

    talk("Run thermal conductivity")
    talk(f"Log will be written to {outfile}")

    sys.stdout = open(outfile, "w")

    rtc.run_thermal_conductivity_in_folder(folder, mesh=q_mesh)
