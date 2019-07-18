""" tools for storing MD trajectories

Logic:
* save md metadata to new trajectory
* append each md step afterwards

"""

import os
import shutil

import numpy as np

from ase import units
from hilde import __version__ as version
from hilde import son
from hilde.fourier import get_timestep
from hilde.helpers.converters import results2dict, dict2atoms, input2dict
from hilde.helpers.converters import dict2json as dumper
from hilde.helpers.hash import hash_atoms
from hilde.helpers import Timer, warn, talk
from hilde.helpers.utils import progressbar


def step2file(atoms, calc=None, file="trajectory.son", append_cell=True, metadata={}):
    """Save the current step

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The structure at the current step
    calc: ase.calculators.calulator.Calculator
        The ASE Calculator for the current run
    file: str or Path
        Path to file to append the current step to
    append_cell: True
        If True add cell to the calculation
    metadata: dict
        the metadata for the calculation, store to atoms.info if possible
    """

    dct = {}
    if metadata:
        for key, val in metadata.items():
            if key in atoms.info and atoms.info[key] == val:
                continue
            else:
                atoms.info.update({"metadata": metadata})
                break

    dct.update(results2dict(atoms, calc, append_cell))

    son.dump(dct, file, dumper=dumper)


def metadata2file(metadata, file="metadata.son"):
    """save metadata to file

    Parameters
    ----------
    metadata: dict
        the metadata to save
    file: str
        filepath to the output file
    """

    if metadata is None:
        metadata = {}

    son.dump({**metadata, "hilde": {"version": version}}, file, is_metadata=True)


def get_hashes_from_trajectory(trajectory, verbose=False):
    """return all hashes from trajectory"""

    try:
        traj = reader(trajectory, verbose=verbose)
    except (FileNotFoundError, KeyError):
        return []

    hashes = []
    for atoms in traj:
        try:
            hashes.append(atoms.info["hash"])
        except (KeyError, AttributeError):
            hashes.append(hash_atoms(atoms))

    return hashes



class Trajectory(list):
    """ A Trajectory is basically a list of Atoms objects with some functionality, e.g.
           - extract and plot several statistics on the MD trajectory
           - convert to other formats like xyz or TDEP """

    def __init__(self, *args, metadata=None):
        """Initializer

        Parameters
        ----------
        metadata: dict
            The metadata for a particular run
        """
        super().__init__(*args)

        if metadata:
            self._metadata = metadata
        else:
            self._metadata = {}

        # a bit of lazy eval
        self._times = None

    @classmethod
    def from_file(cls, file):
        """ Read trajectory from file """
        trajectory = reader(file)
        return trajectory

    @property
    def metadata(self):
        """ Return metadata """
        return self._metadata

    @metadata.setter
    def metadata(self, metadata):
        """Set the metadata"""
        assert isinstance(metadata, dict)
        self._metadata = metadata

    #     fkdev: Might be useful?
    #     @property
    #     def ref_atoms(self):
    #         """ Reference atoms object for computing displacements etc """
    #         if "supercell" in self.metadata:
    #             return dict2atoms(self.metadata["supercell"]["atoms"])
    #         else:
    #             return self[0]

    @property
    def primitive(self):
        """ Return the primitive cell if it is there """
        if "primitive" in self.metadata:
            return dict2atoms(self.metadata["primitive"]["atoms"])
        warn("primitive cell not provided in trajectory metadata")

    @primitive.setter
    def primitive(self, atoms):
        """ Set the supercell atoms object """
        dct = input2dict(atoms)

        self.metadata["primitive"] = dct
        talk(".. primitive added to metadata.")

    @property
    def supercell(self):
        """ Return the supercell if it is there """
        if "supercell" in self.metadata:
            return dict2atoms(self.metadata["supercell"]["atoms"])
        warn("supercell not provided in trajectory metadata")

    @supercell.setter
    def supercell(self, atoms):
        """ Set the supercell atoms object """
        dct = input2dict(atoms)

        self.metadata["supercell"] = dct
        talk(".. supercell added to metadata.")

    @property
    def times(self):
        """ return the times as numpy array in fs"""
        if self._times is None:
            try:
                fs = self.metadata["MD"]["fs"]
            except KeyError:
                warn("time unit not found in trajectory metadata, use ase.units.fs")
                fs = units.fs

            times = np.array([a.info["nsteps"] * a.info["dt"] / fs for a in self])
            self._times = times

        return self._times

    @property
    def timestep(self):
        """ return the timestep in fs"""
        return get_timestep(self.times)

    @property
    def temperatures(self):
        """ return the temperatues as 1d array """
        return np.array([a.get_temperature() for a in self])

    @property
    def velocities(self):
        """return the velocities as [N_t, N_a, 3] array"""
        return np.array([a.get_velocities() for a in self])

    @property
    def stress(self):
        """retunr the stress as [N_t, N_a, 3, 3] array"""
        zeros = np.zeros([3, 3])
        stresses = []
        for a in self:
            if "stress" in a.calc.results:
                stresses.append(a.get_stress(voigt=False))
            else:
                stresses.append(zeros)
        return np.array(stresses)

    def get_pressure(self, GPa=False):
        """return the pressure as [N_t] array"""
        pressure = np.array([-1 / 3 * np.trace(stress) for stress in self.stress])
        if GPa:
            pressure /= units.GPa
        return pressure

    @property
    def pressure(self):
        """return the pressure as [N_t] array"""
        return self.get_pressure()

    def with_result(self, result="stresses"):
        """return new trajectory with atoms object that have specific result computed"""
        atoms_w_result = [a for a in self if result in a.calc.results]
        new_traj = Trajectory(atoms_w_result, metadata=self.metadata)
        return new_traj

    @property
    def with_stresses(self):
        """return new trajectory with atoms that have stresses computed"""
        return self.with_result("stresses")

    def clean_drift(self):
        """ Clean constant drift CAUTION: respect ASE time unit correctly! """

        timer = Timer("Clean trajectory from constant drift")

        p_drift = np.mean([a.get_momenta().sum(axis=0) for a in self], axis=0)

        talk(f".. drift momentum is {p_drift}")

        for atoms, time in zip(self, self.times):
            atoms.set_momenta(atoms.get_momenta() - p_drift / len(atoms))

            # the displacement
            disp = p_drift / atoms.get_masses().sum() * time
            atoms.positions = atoms.positions - disp

        timer("velocities and positions cleaned from drift")

    def write(self, file="trajectory.son"):
        """Write to son file

        Parameters
        ----------
        file: str
            path to trajecotry son file
        """

        timer = Timer(f"Write trajectory to {file}")

        temp_file = "temp.son"

        # check for file and make backup
        if os.path.exists(file):
            ofile = f"{file}.bak"
            shutil.copy(file, ofile)
            talk(f".. {file} copied to {ofile}")

        metadata2file(self.metadata, temp_file)

        prefix = f"Write to {temp_file}:"
        for elem in progressbar(self, prefix=prefix):
            son.dump(results2dict(elem), temp_file)

        shutil.move(temp_file, file)

        timer()

    def to_xyz(self, file="positions.xyz"):
        """Write positions to simple xyz file for e.g. viewing with VMD

        Parameters
        ----------
        file: str
            path to trajecotry xyz file
        """
        from ase.io.xyz import simple_write_xyz

        with open(file, "w") as fo:
            simple_write_xyz(fo, self)

    def to_tdep(self, folder=".", skip=1):
        """Convert to TDEP infiles for direct processing

        Parameters
        ----------
        folder: str or Path
            Directory to store tdep files
        skip: int
            Number of structures to skip
        """
        from pathlib import Path
        from contextlib import ExitStack

        folder = Path(folder)
        folder.mkdir(exist_ok=True)

        talk(f"Write tdep input files to {folder}:")

        # meta
        n_atoms = len(self[0])
        n_steps = len(self) - skip
        try:
            dt = self.metadata["MD"]["timestep"] / self.metadata["MD"]["fs"]
            T0 = self.metadata["MD"]["temperature"] / units.kB
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
        if self.primitive:
            fname = folder / "infile.ucposcar"
            self.primitive.write(str(fname), **write_settings)
            talk(f".. {fname} written.")
        if self.supercell:
            fname = folder / "infile.ssposcar"
            self.supercell.write(str(fname), **write_settings)
            talk(f".. {fname} written.")

        with ExitStack() as stack:
            pdir = folder / "infile.positions"
            fdir = folder / "infile.forces"
            sdir = folder / "infile.stat"
            fp = stack.enter_context(pdir.open("w"))
            ff = stack.enter_context(fdir.open("w"))
            fs = stack.enter_context(sdir.open("w"))

            for ii, atoms in enumerate(self[skip:]):
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

    def get_average_displacements(self, ref_atoms=None, window=-1):
        """Return averaged displacements

        Parameters
        ----------
        ref_atoms: ase.atoms.Atoms
            reference structure for undisplaced system
        window: int
            This does nothing? I think it is supposed to define which steps to start/end the analysis on

        Returns
        -------
        avg_displacement: np.ndarray
            The average displacements of all the atoms in self
        """

        from hilde.harmonic_analysis.displacements import get_dR

        # reference atoms
        if not ref_atoms:
            if "supercell" in self.metadata:
                ref_atoms = dict2atoms(self.metadata["supercell"]["atoms"])
            else:
                ref_atoms = self[0]

        # this will hold the averaged displacement
        avg_displacement = np.zeros_like(ref_atoms.get_positions())

        weigth = 1 / len(self)

        for atoms in self:
            avg_displacement += weigth * get_dR(ref_atoms, atoms)

        return avg_displacement

    def get_average_positions(self, ref_atoms=None, window=-1, wrap=False):
        """ Return averaged positions

        Parameters
        ----------
        ref_atoms: ase.atoms.Atoms
            reference structure for undisplaced system
        window: int
            This does nothing? I think it is supposed to define which steps to start/end the analysis on
        wrap: bool
            If True wrap all the atoms to be within the unit cell

        Returns
        -------
        np.ndarray
            The average positions of all the atoms in self
        """
        # reference atoms
        if not ref_atoms:
            if "supercell" in self.metadata:
                ref_atoms = dict2atoms(self.metadata["supercell"]["atoms"])
            else:
                ref_atoms = self[0]

        avg_displacement = self.get_average_displacements(
            ref_atoms=ref_atoms, window=window
        )

        avg_atoms = ref_atoms.copy()
        avg_atoms.positions += avg_displacement

        if wrap:
            avg_atoms.wrap()

        return avg_atoms.get_positions()

    def get_hashes(self, verbose=False):
        """return all hashes from trajectory"""

        hashes = []
        for atoms in self:
            try:
                hashes.append(atoms.info["hash"])
            except (KeyError, AttributeError):
                hashes.append(hash_atoms(atoms))

        return hashes
