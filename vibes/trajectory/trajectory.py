"""the vibes.Trajectory class"""

import numpy as np

from ase import units, Atoms
from ase.calculators.calculator import PropertyNotImplementedError
from vibes.fourier import get_timestep
from vibes.helpers.converters import dict2atoms, atoms2dict
from vibes.helpers.hash import hash_atoms
from vibes.helpers import warn, lazy_property
from vibes.helpers.utils import progressbar
from vibes.helpers.displacements import get_dR
from vibes.anharmonicity_score import get_sigma
from . import (
    io,
    heat_flux as hf,
    talk,
    Timer,
    dataset as xr,
    analysis as al,
    _prefix,
    _fc_key,
)


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
        self._reference_atoms = None
        self._average_atoms = None
        self._heat_flux = None
        self._avg_heat_flux = None
        self._volume = None
        self._forces_harmonic = None
        self._displacements = None
        self._force_constants_remapped = None
        if _fc_key not in self._metadata:
            self._metadata[_fc_key] = None

    @classmethod
    def read(cls, file="trajectory.son", **kwargs):
        """ Read trajectory from file """
        trajectory = io.reader(file, **kwargs)
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
    def reference_atoms(self):
        """Reference atoms object for computing displacements etc"""
        if not self._reference_atoms:
            if self.supercell:
                self._reference_atoms = self.supercell.copy()
            else:
                warn("No supercell found, return first Atoms in trajectory", level=1)
                self._reference_atoms = self[0].copy()
        return self._reference_atoms

    @reference_atoms.setter
    def reference_atoms(self, atoms):
        assert isinstance(atoms, Atoms)
        self._reference_atoms = atoms

    # ref_atoms legacy
    @property
    def ref_atoms(self):
        return self.reference_atoms

    @ref_atoms.setter
    def ref_atoms(self, atoms):
        self.reference_atoms = atoms

    @property
    def primitive(self):
        """ Return the primitive atoms if it is there """
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
        return self.reference_atoms.get_chemical_symbols()

    @property
    def masses(self):
        """return masses in AMU"""
        return self.reference_atoms.get_masses()

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

        try:
            return np.array([a.info["nsteps"] * a.info["dt"] / fs for a in self])
        except KeyError:
            warn("no time steps found, return time as index", level=1)
            return np.arange(len(self))

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
        return self.reference_atoms.get_positions()

    @lazy_property
    def positions(self):
        """return the positions as [N_t, N_a, 3] array"""
        return np.array([a.get_positions() for a in self])

    @lazy_property
    def velocities(self):
        """return the velocities as [N_t, N_a, 3] array"""
        return np.array([a.get_velocities() for a in self])

    @lazy_property
    def momenta(self):
        """return the velocities as [N_t, N_a, 3] array"""
        return np.array([a.get_momenta() for a in self])

    @lazy_property
    def forces(self):
        """return the forces as [N_t, N_a, 3] array"""
        return np.array([a.get_forces() for a in self])

    @property
    def force_constants(self):
        """return (reduced) force constants or warn if not set"""
        fc = self.metadata[_fc_key]
        if any(x is None for x in (fc, self.primitive, self.supercell)):
            warn("`trajectory.force_constants` not set, return None")
        else:
            fc = np.asarray(fc)
            Np, Na = len(self.primitive), len(self.supercell)
            assert fc.shape == (Np, Na, 3, 3), fc.shape

            return fc

    def set_force_constants(self, fc):
        """set force constants"""
        Np, Na = len(self.primitive), len(self.supercell)
        # assert fc.shape == 2 * (3 * len(self.supercell),), fc.shape
        assert fc.shape == (Np, Na, 3, 3), fc.shape
        self.metadata[_fc_key] = fc

    @property
    def force_constants_remapped(self):
        """return remapped force constants [3 * Na, 3 * Na]"""
        if self._force_constants_remapped is None:
            fc = self.force_constants
            if fc is not None:
                uc, sc = self.primitive, self.supercell
                from vibes.phonopy.utils import remap_force_constants

                fc = remap_force_constants(fc, uc, sc, two_dim=True, symmetrize=True)

                self._force_constants_remapped = fc
        return self._force_constants_remapped

    def set_force_constants_remapped(self, fc):
        """set remapped force constants"""
        Na = len(self.supercell)
        assert fc.shape == (3 * Na, 3 * Na), fc.shape
        self._force_constants_remapped = fc

    def set_forces_harmonic(self, force_constants=None, average_reference=False):
        """return harmonic force computed from force_constants"""
        if average_reference:
            talk("Compute harmonic force constants with average positions as reference")
            self.reference_atoms = self.average_atoms

        if force_constants is None:
            force_constants = self.force_constants_remapped

        timer = Timer("Set harmonic forces")

        forces_ha = [-force_constants @ d.flatten() for d in self.displacements]
        timer()

        self._forces_harmonic = np.array(forces_ha).reshape(self.positions.shape)

    @property
    def forces_harmonic(self):
        """return harmonic forces, None if not set via `set_force_constants`"""
        if self._forces_harmonic is None:
            if self.force_constants_remapped is not None:
                self.set_forces_harmonic()

        return self._forces_harmonic

    @lazy_property
    def kinetic_energy(self):
        """return the kinetic energy as [N_t] array"""
        return np.array([a.get_kinetic_energy() for a in self])

    @lazy_property
    def potential_energy(self):
        """return the potential energy as [N_t] array"""
        return np.array([a.get_potential_energy() for a in self])

    @lazy_property
    def potential_energy_harmonic(self):
        return -0.5 * (self.forces_harmonic * self.displacements).sum(axis=(1, 2))

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

        return atomic_stresses

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
    def sigma(self):
        """return sigma_A"""
        if self.forces_harmonic is not None:
            x = self.forces
            y = self.forces_harmonic
            return get_sigma(x, y)

    @property
    def sigma_per_sample(self):
        """return sigma_A per time step"""
        if self.forces_harmonic is not None:
            x = self.forces
            y = self.forces_harmonic
            return get_sigma(x, y, axis=(1, 2))

    @lazy_property
    def dataset(self):
        """return data as xarray.Dataset

        Contains:
            positions, velocities, forces, stress, pressure, temperature
        """
        return xr.get_trajectory_dataset(self)

    def get_heat_flux_data(self, only_flux=False):
        """return heat flux and other data as xarray.Dataset

        Contains:
            heat_flux, avg_heat_flux, velocities, forces, pressure, temperature
        Metadata:
            volume, symbols, masses, flattend reference positions
        """
        return xr.get_heat_flux_dataset(self, only_flux=only_flux)

    @property
    def heat_flux_dataset(self):
        """== self.get_heat_flux_data()"""
        return self.get_heat_flux_data()

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
        """Write to son or nc file

        Args:
            file: path to trajecotry son or nc file
        """
        io.write(self, file=file)

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

    def to_db(self, database):
        """Convert to ase database

        Args:
            database: Filename or address of database

        """
        io.to_db(self, database)

    def set_displacements(self):
        """calculate the displacements for `reference_atoms`"""
        if not self.supercell:
            # warn("Supercell not set, let us stop here.", level=2)
            warn("SUPERCELL NOT SET, compute w.r.t to reference atoms", level=1)

        atoms_ideal = self.reference_atoms

        # use get_dR
        list1 = []
        msg = f"[{_prefix}]   compute displacements: "
        for atoms in progressbar(self, prefix=msg):
            list1.append(get_dR(atoms, atoms_ideal))
        list1 = np.asarray(list1)

        self.displacements = list1

    @property
    def displacements(self):
        """cached version of `get_displacements`"""
        if self._displacements is None:
            self.set_displacements()
        return self._displacements

    @displacements.setter
    def displacements(self, displacement_array):
        assert np.shape(displacement_array) == np.shape(self.positions)
        self._displacements = displacement_array

    def get_average_displacements(self, window=-1):
        """Return averaged displacements

        Args:
            window: This does nothing
        Returns:
            array: The average displacements of all the atoms in self
        """

        displacements = self.displacements

        weight = 1  # 1 / len(displacements)

        # this will hold the averaged displacement
        avg_displacement = weight * displacements.mean(axis=0)

        return avg_displacement

    def get_average_positions(self, window=-1, wrap=False):
        """ Return averaged positions

        Args:
            window: This does nothing
            wrap: If True wrap all the atoms to be within the unit cell

        Returns:
            np.ndarray: The average positions of all the atoms in self
        """
        # reference atoms
        ref_atoms = self.reference_atoms

        avg_displacement = self.get_average_displacements(window=window)

        avg_atoms = ref_atoms.copy()
        avg_atoms.positions += avg_displacement

        if wrap:
            avg_atoms.wrap()

        return avg_atoms.get_positions()

    @property
    def average_atoms(self):
        """Atoms object with averaged positions"""
        if not self._average_atoms:
            self.set_average_reference()
        return self._average_atoms

    def set_average_reference(self):
        talk(f"(Re-)set average positions")
        avg_atoms = self.reference_atoms.copy()
        avg_atoms.positions = self.get_average_positions()
        self._average_atoms = avg_atoms

    def get_hashes(self, verbose=False):
        """return all hashes from trajectory"""

        hashes = []
        for atoms in self:
            hashes.append(hash_atoms(atoms))

        return hashes

    def summarize(self, vebose=False):
        """give a summary of relevant statistics"""

        DS = self.dataset

        al.pressure(DS.pressure)
