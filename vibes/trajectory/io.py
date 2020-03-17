"""Trajectory File I/O"""

import json
import os
import shutil
from pathlib import Path

import numpy as np
import xarray as xr
from ase import Atoms, units
from ase.calculators.calculator import PropertyNotImplementedError
from ase.calculators.singlepoint import SinglePointCalculator
from vibes import __version__ as version
from vibes import io, keys, son
from vibes.helpers import warn
from vibes.helpers.converters import dict2atoms, dict2json, results2dict
from vibes.helpers.utils import progressbar

from .utils import Timer, talk


def step2file(atoms, calc=None, file="trajectory.son", metadata={}):
    """Save the current step

    Args:
        atoms: The structure at the current step
        calc: The ASE Calculator for the current run
        file: Path to file to append the current step to
        metadata: the metadata for the calculation, store to atoms.info if possible
    """

    dct = {}
    if metadata:
        for key, val in metadata.items():
            if key in atoms.info and atoms.info[key] == val:
                continue
            if key not in atoms.info:
                atoms.info[key] = val
            else:
                atoms.info.update({"metadata": metadata})
                break

    dct.update(results2dict(atoms, calc))

    son.dump(dct, file, dumper=dict2json)


def metadata2file(metadata, file="metadata.son"):
    """save metadata to file

    Args:
        metadata: the metadata to save
        file: filepath to the output file
    """

    if metadata is None:
        metadata = {}

    son.dump({**metadata, "vibes": {"version": version}}, file, is_metadata=True)


def write(trajectory, file="trajectory.son"):
    """Write to son file

    Args:
        file: path to trajecotry son or netcdf file
    """
    from .dataset import get_trajectory_dataset

    timer = Timer(f"Write trajectory to {file}")

    if Path(file).suffix == ".nc":
        dataset = get_trajectory_dataset(trajectory, metadata=True)
        dataset.to_netcdf(file)
        timer()
        return True

    temp_file = "temp.son"

    # check for file and make backup
    if os.path.exists(file):
        ofile = f"{file}.bak"
        shutil.copy(file, ofile)
        talk(f".. {file} copied to {ofile}")

    metadata2file(trajectory.metadata, temp_file)

    talk(f"Write to {temp_file}")
    for elem in progressbar(trajectory):
        son.dump(results2dict(elem), temp_file)

    shutil.move(temp_file, file)

    timer()


def reader(
    file="trajectory.son",
    get_metadata=False,
    fc_file=None,
    with_stresses=False,
    verbose=True,
    single_point_calc=True,
):
    """Convert information in file to Trajectory

    Args:
        file: Trajectory file to read the structures from
        get_metadata: If True return the metadata
        fc_file: force constants file
        with_stresses: Return only the atoms with stresses computed
        verbose: If True print more information to the screen

    Returns:
        trajectory: The trajectory from the file
        metadata: The metadata for the trajectory
    """
    from vibes.trajectory.trajectory import Trajectory

    timer = Timer(f"Parse `{file}`")

    if Path(file).suffix == ".nc":
        trajectory = read_netcdf(file)
        timer()

        return trajectory

    try:
        metadata, pre_trajectory = son.load(file, verbose=verbose)
    except json.decoder.JSONDecodeError:
        metadata, pre_trajectory = son.load(file, verbose=verbose)

    # legacy of trajectory.yaml
    if metadata is None:
        msg = f"metadata in {file} appears to be empty, assume old convention w/o === "
        msg += f"was used. Let's see"
        warn(msg, level=1)
        metadata = pre_trajectory.pop(0)

    pre_calc_dict = metadata["calculator"]
    pre_atoms_dict = metadata["atoms"]

    if "numbers" in pre_atoms_dict and "symbols" in pre_atoms_dict:
        del pre_atoms_dict["symbols"]

    if "MD" in metadata:
        md_metadata = metadata["MD"]

    if not pre_trajectory:
        if get_metadata:
            talk(".. trajectory empty, return (Trajectory([]), metadata)")
            return Trajectory([]), metadata
        talk(".. trajectory empty, return Trajectory([])")
        return Trajectory([])

    trajectory = Trajectory(metadata=metadata)
    talk(".. create atoms")
    for obj in progressbar(pre_trajectory):

        atoms_dict = {**pre_atoms_dict, **obj["atoms"]}

        # remember that the results need to go to a dedicated results dict in calc
        calc_dict = {**pre_calc_dict, "results": obj["calculator"]}

        atoms = dict2atoms(atoms_dict, calc_dict, single_point_calc)

        # info
        if "MD" in metadata:
            if "dt" in atoms.info:
                atoms.info["dt_fs"] = atoms.info["dt"] / md_metadata["fs"]
        elif "info" in obj:
            info = obj["info"]
            atoms.info.update(info)

        # compatibility with older trajectories
        if "MD" in obj:
            atoms.info.update(obj["MD"])

        # preserve metadata
        if "metadata" in obj:
            atoms.info.update({"metadata": obj["metadata"]})

        trajectory.append(atoms)

    timer("done")

    if with_stresses:
        talk(".. return only atoms with `stresses` computed")
        trajectory = trajectory.with_stresses

    if fc_file:
        fc = io.parse_force_constants(fc_file, two_dim=False)
        trajectory.set_force_constants(fc)

    if get_metadata:
        return trajectory, metadata
    return trajectory


def to_tdep(trajectory, folder=".", skip=1):
    """Convert to TDEP infiles for direct processing

    Args:
        folder: Directory to store tdep files
        skip: Number of structures to skip
    """
    from pathlib import Path
    from contextlib import ExitStack

    folder = Path(folder)
    folder.mkdir(exist_ok=True)

    talk(f"Write tdep input files to {folder}:")

    # meta
    n_atoms = len(trajectory[0])
    n_steps = len(trajectory) - skip
    try:
        dt = trajectory.metadata["MD"]["timestep"] / trajectory.metadata["MD"]["fs"]
        T0 = trajectory.metadata["MD"]["temperature"] / units.kB
    except KeyError:
        dt = 1.0
        T0 = 0

    lines = [f"{n_atoms}", f"{n_steps}", f"{dt}", f"{T0}"]

    fname = folder / "infile.meta"

    with fname.open("w") as fo:
        fo.write("\n".join(lines))
        talk(f".. {fname} written.")

    # supercell and fake unit cell
    write_settings = {"format": "vasp", "direct": True, "vasp5": True}
    if trajectory.primitive:
        fname = folder / "infile.ucposcar"
        trajectory.primitive.write(str(fname), **write_settings)
        talk(f".. {fname} written.")
    if trajectory.supercell:
        fname = folder / "infile.ssposcar"
        trajectory.supercell.write(str(fname), **write_settings)
        talk(f".. {fname} written.")

    with ExitStack() as stack:
        pdir = folder / "infile.positions"
        fdir = folder / "infile.forces"
        sdir = folder / "infile.stat"
        fp = stack.enter_context(pdir.open("w"))
        ff = stack.enter_context(fdir.open("w"))
        fs = stack.enter_context(sdir.open("w"))

        for ii, atoms in enumerate(trajectory[skip:]):
            # stress and pressure in GPa
            try:
                stress = atoms.get_stress(voigt=True) / units.GPa
                pressure = -1 / 3 * sum(stress[:3])
            except PropertyNotImplementedError:
                stress = np.zeros(6)
                pressure = 0.0
            e_tot = atoms.get_total_energy()
            e_kin = atoms.get_kinetic_energy()
            e_pot = e_tot - e_kin
            temp = atoms.get_temperature()

            for spos in atoms.get_scaled_positions():
                fp.write("{} {} {}\n".format(*spos))

            for force in atoms.get_forces():
                ff.write("{} {} {}\n".format(*force))

            stat = (
                f"{ii:5d} {ii*dt:10.2f} {e_tot:20.8f} {e_pot:20.8f} "
                f"{e_kin:20.15f} {temp:20.15f} {pressure:20.15f} "
            )
            stat += " ".join([str(s) for s in stress])

            fs.write(f"{stat}\n")

    talk(f".. {sdir} written.")
    talk(f".. {pdir} written.")
    talk(f".. {fdir} written.")


def to_db(trajectory, database):
    """Write vibes trajectory as ase database

    Always creates a new database. Database type
    is always inferred from the filename. Metadata
    is carried over to the ase database.

    Please be aware that ase.db indices start from
    1, not from 0 as usual.

    Args:
        trajectory: Trajectory instance
        database: Filename or address of database

    """
    from ase.db import connect

    timer = Timer(f"Save as ase database {database}")

    with connect(database, append=False) as db:
        talk("write db")
        for atoms in progressbar(trajectory):
            db.write(atoms, data={"info": atoms.info})

    # metadata can only be written *after* the database exists
    with connect(database) as db:
        db.metadata = trajectory.metadata

    timer("done")


def read_netcdf(file="trajectory.nc"):
    """read `trajectory.nc` and return Trajectory"""
    from vibes.trajectory.trajectory import Trajectory

    DS = xr.open_dataset(file)
    attrs = DS.attrs

    # check mandatory keys
    assert keys.reference_atoms in attrs
    assert "positions" in DS
    assert "velocities" in DS
    assert "forces" in DS
    assert keys.energy_potential in DS

    atoms_dict = json.loads(attrs[keys.reference_atoms])

    # metadata
    metadata = None
    if keys.metadata in attrs:
        metadata = json.loads(attrs[keys.metadata])

    # popping `velocities` is obsolete if
    # https://gitlab.com/ase/ase/merge_requests/1563
    # is accepted
    velocities = atoms_dict.pop("velocities", None)

    ref_atoms = Atoms(**atoms_dict)

    if velocities is not None:
        ref_atoms.set_velocities(velocities)

    positions = DS.positions.data
    velocities = DS.velocities.data
    forces = DS.forces.data
    potential_energy = DS[keys.energy_potential].data

    if "stress" in DS:
        stress = DS.stress.data
    else:
        stress = len(positions) * [None]

    traj = []
    for p, v, f, e, s in zip(positions, velocities, forces, potential_energy, stress):
        atoms = ref_atoms.copy()
        atoms.set_positions(p)
        atoms.set_velocities(v)

        results = {"energy": e, "forces": f}
        if s is not None:
            results.update({"stress": s})

        calc = SinglePointCalculator(atoms, **results)
        atoms.calc = calc

        traj.append(atoms)

    trajectory = Trajectory(traj, metadata=metadata)
    trajectory.displacements = DS.displacements.data
    trajectory.times = DS.time.data

    return trajectory
