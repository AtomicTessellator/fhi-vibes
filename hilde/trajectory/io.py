"""Trajectory File I/O"""

import json
import numpy as np
from ase import units
from hilde import son
from hilde import __version__ as version
from hilde.helpers.converters import dict2atoms
from hilde.helpers import Timer, warn, talk
from hilde.helpers.converters import results2dict
from hilde.helpers.converters import dict2json as dumper
from hilde.helpers.utils import progressbar
from hilde.trajectory.trajectory import Trajectory


def step2file(atoms, calc=None, file="trajectory.son", append_cell=True, metadata={}):
    """Save the current step

    Args:
        atoms: The structure at the current step
        calc: The ASE Calculator for the current run
        file: Path to file to append the current step to
        append_cell: If True add cell to the calculation
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

    dct.update(results2dict(atoms, calc, append_cell))

    son.dump(dct, file, dumper=dumper)


def metadata2file(metadata, file="metadata.son"):
    """save metadata to file

    Args:
        metadata: the metadata to save
        file: filepath to the output file
    """

    if metadata is None:
        metadata = {}

    son.dump({**metadata, "hilde": {"version": version}}, file, is_metadata=True)


def reader(file="trajectory.son", get_metadata=False, verbose=True):
    """Convert information in file to Trajectory

    Args:
        trajectory: Trajectory file to pull the structures from
        get_metadata: If True return the metadata
        verbose: If True print more information to the screen

    Returns:
        trajectory: The trajectory from the file
        metadata: The metadata for the trajectory
    """
    timer = Timer(f"Parse trajectory")

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
            talk(".. trajectory empty, return ([], metadata)")
            return [], metadata
        talk(".. trajectory empty, return []")
        return []

    trajectory = Trajectory(metadata=metadata)
    prefix = ".. create atoms: "
    for obj in progressbar(pre_trajectory, prefix=prefix):

        atoms_dict = {**pre_atoms_dict, **obj["atoms"]}

        # remember that the results need to go to a dedicated results dict in calc
        calc_dict = {**pre_calc_dict, "results": obj["calculator"]}

        atoms = dict2atoms(atoms_dict, calc_dict)

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
            except:
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
