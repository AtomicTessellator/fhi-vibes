"""the hilde.Trajectory class"""

import os
import shutil

import numpy as np

from ase import units, Atoms
from ase.calculators.calculator import PropertyNotImplementedError
from hilde import son
from hilde.fourier import get_timestep
from hilde.helpers.converters import results2dict, dict2atoms, atoms2dict
from hilde.helpers.hash import hash_atoms
from hilde.helpers import warn, lazy_property
from hilde.helpers.utils import progressbar
from . import io, heat_flux as hf, talk, Timer, dataarray as xr, analysis as al


class Trajectory(list):
    """ A Trajectory is basically a list of Atoms objects with some functionality, e.g.
           - extract and plot several statistics on the MD trajectory
           - convert to other formats like xyz or TDEP """

    def __init__(self, *args, metadata=None):
        """Initializer

        Args:
            metadata: The metadata for a particular run
        """
        super().__init__(*args)

        if metadata:
            self._metadata = metadata
        else:
            self._metadata = {}

        # lazy eval where @lazy_eval is not applicable
        self._supercell = None
        self._heat_flux = None
        self._avg_heat_flux = None
        self._volume = None

    @classmethod
    def from_file(cls, file):
        """ Read trajectory from file """
        trajectory = io.reader(file)
        return trajectory

    def __getitem__(self, key):
        """returns `trajectory[key]` as Atoms object or new Trajectory instance"""
        temp = super().__getitem__(key)
        if isinstance(temp, Atoms):
            return temp
        metadata = self.metadata
        new_trajectory = Trajectory(temp, metadata=metadata)
        return new_trajectory

    @property
    def metadata(self):
        """ Return metadata """
        return self._metadata

    @metadata.setter
    def metadata(self, metadata):
        """Set the metadata"""
        assert isinstance(metadata, dict)
        self._metadata = metadata

    @property
    def ref_atoms(self):
        """Reference atoms object for computing displacements etc"""
        if self.supercell:
            return self.supercell
        warn("No supercell found, return first Atoms in trajectory", level=1)
        return self[0]

    @property
    def primitive(self):
        """ Return the primitive cell if it is there """
        if "primitive" in self.metadata:
            dct = self.metadata["primitive"]
            if "atoms" in dct:
                dct = dct["atoms"]
            return dict2atoms(dct)
        warn("primitive cell not provided in trajectory metadata")

    @primitive.setter
    def primitive(self, atoms):
        """ Set the supercell atoms object """
        assert isinstance(atoms, Atoms)
        dct = atoms2dict(atoms)
        self.metadata["primitive"] = dct
        talk(".. primitive added to metadata.")

    @property
    def supercell(self):
        """ Return the supercell if it is there """
        if not self._supercell:
            if "supercell" in self.metadata:
                dct = self.metadata["supercell"]
                if "atoms" in dct:
                    dct = dct["atoms"]
                self._supercell = dict2atoms(dct)
            else:
                warn("supercell not provided in trajectory metadata")
        return self._supercell

    @supercell.setter
    def supercell(self, atoms):
        """ Set the supercell atoms object """
        assert isinstance(atoms, Atoms)
        dct = atoms2dict(atoms)
        self.metadata["supercell"] = dct
        talk(".. supercell added to metadata.")
        # also add as attribute
        self._supercell = atoms

    @property
    def symbols(self):
        """return chemical symbols"""
        return self.ref_atoms.get_chemical_symbols()

    @property
    def masses(self):
        """return masses in AMU"""
        return self.ref_atoms.get_masses()

    @property
    def masses_dict(self):
        """return masses in AMU as dictionary"""
        masses_dict = {}
        for sym, mass in zip(self.symbols, self.masses):
            masses_dict[sym] = mass
        return masses_dict

    @property
    def volume(self):
        """return averaged volume"""
        if not self._volume:
            volumes = [a.get_volume() for a in self]
            self._volume = np.mean(volumes).squeeze()
        return self._volume

    @lazy_property
    def times(self):
        """ return the times as numpy array in fs"""
        try:
            fs = self.metadata["MD"]["fs"]
        except KeyError:
            warn("time unit not found in trajectory metadata, use ase.units.fs")
            fs = units.fs

        return np.array([a.info["nsteps"] * a.info["dt"] / fs for a in self])

    @property
    def timestep(self):
        """ return the timestep in fs"""
        return get_timestep(self.times)

    @lazy_property
    def temperatures(self):
        """return the temperatues as 1d array"""
        return np.array([a.get_temperature() for a in self])

    @property
    def ref_positions(self):
        """return reference positions"""
        return self.ref_atoms.get_positions()

    @lazy_property
    def positions(self):
        """return the positions as [N_t, N_a, 3] array"""
        return np.array([a.get_positions() for a in self])

    @lazy_property
    def velocities(self):
        """return the velocities as [N_t, N_a, 3] array"""
        return np.array([a.get_velocities() for a in self])

    @lazy_property
    def forces(self):
        """return the forces as [N_t, N_a, 3] array"""
        return np.array([a.get_forces() for a in self])

    @lazy_property
    def kinetic_energy(self):
        """return the kinetic energy as [N_t] array"""
        return np.array([a.get_kinetic_energy() for a in self])

    @lazy_property
    def potential_energy(self):
        """return the potential energy as [N_t] array"""
        return np.array([a.get_potential_energy() for a in self])

    @lazy_property
    def stress(self):
        """return the stress as [N_t, 3, 3] array"""
        zeros = np.zeros([3, 3])
        stresses = []
        for a in self:
            if "stress" in a.calc.results:
                stresses.append(a.get_stress(voigt=False))
            else:
                stresses.append(zeros)

        return np.array(stresses)

    @lazy_property
    def stresses(self):
        """return the atomic stress as [N_t, N_a, 3, 3] array"""
        atomic_stresses = []
        for a in self:
            V = a.get_volume()
            try:
                atomic_stress = a.get_stresses() / V
            except PropertyNotImplementedError:
                msg = "Trajectory contains atoms without stresses computed. "
                msg += "Use `trajectory.with_stresses`"
                warn(msg, level=2)
            atomic_stresses.append(atomic_stress)

        return atomic_stress

    @property
    def heat_flux(self):
        """return heat flux as [N_t, N_a, 3] array"""
        if self._heat_flux is None:
            self._heat_flux, self._avg_heat_flux = hf.get_heat_flux(self)
        return self._heat_flux

    @property
    def avg_heat_flux(self):
        """return heat flux FROM AVERAGED STRESS as [N_t, N_a, 3] array"""
        if self._avg_heat_flux is None:
            self._heat_flux, self._avg_heat_flux = hf.get_heat_flux(self)
        return self._avg_heat_flux

    def get_pressure(self, GPa=False):
        """return the pressure as [N_t] array"""
        pressure = np.array([-1 / 3 * np.trace(stress) for stress in self.stress])
        assert len(pressure) == len(self)

        if GPa:
            pressure /= units.GPa
        return pressure

    @lazy_property
    def pressure(self):
        """return the pressure as [N_t] array"""
        return self.get_pressure()

    @property
    def dataset(self):
        """return data as xarray.Dataset

        Contains:
            positions, velocities, forces, stress, pressure, temperature
        """
        return xr.get_trajectory_data(self)

    @property
    def heat_flux_dataset(self):
        """return heat flux and other data as xarray.Dataset

        Contains:
            heat_flux, avg_heat_flux, velocities, forces, pressure, temperature
        Metadata:
            volume, symbols, masses, flattend reference positions
        """
        return xr.get_heat_flux_data(self)

    def with_result(self, result="stresses"):
        """return new trajectory with atoms object that have specific result computed"""
        atoms_w_result = [a for a in self if result in a.calc.results]
        new_traj = Trajectory(atoms_w_result, metadata=self.metadata)
        return new_traj

    @property
    def with_stresses(self):
        """return new trajectory with atoms that have stresses computed"""
        return self.with_result("stresses")

    def discard(self, first=0, last=0):
        """discard atoms before FIRST and after LAST and return as new Trajectory"""
        n = len(self)
        part = self[first : n - last]
        talk(f"Discard first {first} atoms")
        talk(f"Discard last  {last} atoms")
        talk(f".. length before: {n}")
        talk(f".. length after:  {len(part)}")
        return Trajectory(part, metadata=self.metadata)

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

        Args:
            file: path to trajecotry son file
        """

        timer = Timer(f"Write trajectory to {file}")

        temp_file = "temp.son"

        # check for file and make backup
        if os.path.exists(file):
            ofile = f"{file}.bak"
            shutil.copy(file, ofile)
            talk(f".. {file} copied to {ofile}")

        io.metadata2file(self.metadata, temp_file)

        prefix = f"Write to {temp_file}:"
        for elem in progressbar(self, prefix=prefix):
            son.dump(results2dict(elem), temp_file)

        shutil.move(temp_file, file)

        timer()

    def to_xyz(self, file="positions.xyz"):
        """Write positions to simple xyz file for e.g. viewing with VMD

        Args:
            file: path to trajecotry xyz file
        """
        from ase.io.xyz import simple_write_xyz

        with open(file, "w") as fo:
            simple_write_xyz(fo, self)

    def to_tdep(self, folder=".", skip=1):
        """Convert to TDEP infiles for direct processing

        Args:
            folder: Directory to store tdep files
            skip: Number of structures to skip
        """
        io.to_tdep(self, folder, skip)

    def get_average_displacements(self, ref_atoms=None, window=-1):
        """Return averaged displacements

        Args:
            ref_atoms: reference structure for undisplaced system
            window: This does nothing
        Returns:
        avg_displacement: The average displacements of all the atoms in self
        """

        from hilde.harmonic_analysis.displacements import get_dR

        # reference atoms
        ref_atoms = self.ref_atoms

        # this will hold the averaged displacement
        avg_displacement = np.zeros_like(ref_atoms.get_positions())

        weigth = 1 / len(self)

        for atoms in self:
            avg_displacement += weigth * get_dR(ref_atoms, atoms)

        return avg_displacement

    def get_average_positions(self, ref_atoms=None, window=-1, wrap=False):
        """ Return averaged positions

        Args:
            ref_atoms: reference structure for undisplaced system
            window: This does nothing
            wrap: If True wrap all the atoms to be within the unit cell

        Returns:
            np.ndarray: The average positions of all the atoms in self
        """
        # reference atoms
        ref_atoms = self.ref_atoms

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

    def summarize(self, vebose=False):
        """give a summary of relevant statistics"""

        DS = self.dataset

        al.pressure(DS.pressure)
