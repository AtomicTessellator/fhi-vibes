""" tools for storing MD trajectories

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import json

from ase import units
from hilde import __version__ as version
from hilde.helpers.converters import results2dict, dict2results
from hilde.helpers.fileformats import to_yaml, from_yaml
from hilde.helpers.hash import hash_atoms


def step2file(atoms, calc, file="trajectory.yaml", append_cell=False):
    """ Save the current step """

    to_yaml(results2dict(atoms, calc, append_cell), file)


def metadata2file(metadata, file="metadata.yaml"):
    """ save metadata to file """

    if metadata is None:
        metadata = {}

    to_yaml({**metadata, "hilde": {"version": version}}, file, mode="w")


def get_hashes_from_trajectory(trajectory):
    """ return all hashes from trajectory """

    try:
        traj = reader(trajectory)
    except (FileNotFoundError, KeyError):
        return []

    hashes = []
    for atoms in traj:
        try:
            hashes.append(atoms.info["hash"])
        except (KeyError, AttributeError):
            hashes.append(hash_atoms(atoms))

    return hashes


def reader(file, get_metadata=False):
    """ convert information in trajectory and metadata files to atoms objects
     and return them """

    try:
        metadata, *pre_trajectory = from_yaml(file, use_json=True)
    except json.decoder.JSONDecodeError:
        metadata, *pre_trajectory = from_yaml(file, use_json=False)

    pre_calc_dict = metadata["calculator"]
    pre_atoms_dict = metadata["atoms"]

    if "numbers" in pre_atoms_dict and "symbols" in pre_atoms_dict:
        del pre_atoms_dict["symbols"]

    if "MD" in metadata:
        md_metadata = metadata["MD"]

    trajectory = Trajectory(metadata=metadata)
    for obj in pre_trajectory:

        atoms_dict = {**pre_atoms_dict, **obj["atoms"]}

        # remember that the results need to go to a dedicated results dict in calc
        calc_dict = {**pre_calc_dict, "results": obj["calculator"]}

        atoms = dict2results(atoms_dict, calc_dict)

        # info
        if "MD" in metadata:
            if "dt" in atoms.info:
                atoms.info["dt"] /= md_metadata["fs"]
        elif "info" in obj:
            info = obj["info"]
            atoms.info.update(info)

        # compatibility with older trajectories
        if "MD" in obj:
            atoms.info.update(obj["MD"])

        trajectory.append(atoms)
    if get_metadata:
        return trajectory, metadata
    return trajectory


class Trajectory(list):
    """ A Trajectory is basically a list of Atoms objects with some functionality, e.g.
           - extract and plot several statistics on the MD trajectory
           - convert to other formats like xyz or TDEP """

    def __init__(self, metadata=None):
        super().__init__()

        self._metadata = metadata

    @property
    def metadata(self):
        """ Return metadata """
        return self._metadata

    def to_xyz(self, file="positions.xyz"):
        """ Write positions to simple xyz file for e.g. viewing with VMD """
        from ase.io.xyz import simple_write_xyz

        with open(file, "w") as fo:
            simple_write_xyz(fo, self)

    def to_tdep(self, folder=".", skip=1):
        """ Convert to TDEP infiles for direct processing """
        from pathlib import Path
        from contextlib import ExitStack

        folder = Path(folder)
        folder.mkdir(exist_ok=True)

        # meta
        n_atoms = len(self[0])
        n_steps = len(self) - skip
        dt = self.metadata["MD"]["timestep"] / self.metadata["MD"]["fs"]
        try:
            T0 = self.metadata["MD"]["temperature"] / units.kB
        except KeyError:
            T0 = 0

        lines = [f"{n_atoms}", f"{n_steps}", f"{dt}", f"{T0}"]
        with (folder / "infile.meta").open("w") as fo:
            fo.write("\n".join(lines))

        # supercell and fake unit cell
        write_settings = {"format": "vasp", "direct": True, "vasp5": True}
        self[0].write(str(folder / "infile.ssposcar"), **write_settings)
        self[0].write(str(folder / "infile.ucposcar"), **write_settings)

        with ExitStack() as stack:
            fp = stack.enter_context((folder / "infile.positions").open("w"))
            ff = stack.enter_context((folder / "infile.forces").open("w"))
            fs = stack.enter_context((folder / "infile.stat").open("w"))

            for ii, atoms in enumerate(self[skip:]):
                # stress and pressure in GPa
                stress = atoms.get_stress(voigt=True) / units.GPa
                pressure = 1 / 3 * sum(stress[:3])
                e_tot = atoms.get_total_energy()
                e_kin = atoms.get_kinetic_energy()
                e_pot = e_tot - e_kin
                temp = atoms.get_temperature()

                for spos in atoms.get_scaled_positions():
                    fp.write("{} {} {}\n".format(*spos))

                for force in atoms.get_forces():
                    ff.write("{} {} {}\n".format(*force))

                stat = f"{ii} {ii*dt} {e_tot} {e_pot} {e_kin} {temp} {pressure} "
                stat += " ".join([str(s) for s in stress])

                fs.write(f"{stat}\n")
