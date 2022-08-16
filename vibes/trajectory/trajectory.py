"""the vibes.Trajectory class"""

import numpy as np
from ase import Atoms, units
from ase.geometry import find_mic
from vibes import keys
from vibes.anharmonicity_score import get_sigma
from vibes.filenames import filenames
from vibes.force_constants import ForceConstants
from vibes.helpers import lazy_property, warn
from vibes.helpers.converters import atoms2dict, dict2atoms
from vibes.helpers.hash import hash_atoms
from vibes.helpers.stress import has_stress
from vibes.helpers.stresses import get_stresses, has_stresses
from vibes.helpers.utils import progressbar
from vibes.helpers.virials import get_virials, has_virials

from .utils import Timer, talk


class Trajectory(list):
    """A Trajectory is basically a list of Atoms objects with some functionality, e.g.
    - extract and plot several statistics on the MD trajectory
    - convert to other formats like xyz or TDEP"""

    def __init__(self, *args, metadata=None, debug=False):
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
        self._times = None
        self._supercell = None
        self._reference_atoms = None
        self._average_atoms = None
        self._heat_flux = None
        self._heat_flux_aux = None
        self._heat_flux_harmonic = None
        self._heat_flux_0_harmonic = None
        self._forces_harmonic = None
        self._displacements = None
        self._force_constants = None
        self._uuid = None
        self._hash_raw = None
        if keys.fc not in self._metadata:
            self._metadata[keys.fc] = None

        self.debug = debug

    @classmethod
    def read(cls, file=filenames.trajectory, **kwargs):
        """Read trajectory from file"""
        from .io import reader

        return reader(file, **kwargs)

    @classmethod
    def from_dataset(cls, dataset):
        """parse from xarray.Dataset"""
        from .io import parse_dataset

        return parse_dataset(dataset)

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
        """Return metadata"""
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
                if self.debug:
                    msg = "No supercell found, return first Atoms in trajectory"
                    warn(msg, level=1)
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

    @property
    def primitive(self):
        """Return the primitive atoms if it is there"""
        if "primitive" in self.metadata:
            dct = self.metadata["primitive"]
            if "atoms" in dct:
                dct = dct["atoms"]
            return dict2atoms(dct)
        if self.debug:
            warn("primitive cell not provided in trajectory metadata")

    @primitive.setter
    def primitive(self, atoms):
        """ Set the primitive atoms object """
        assert isinstance(atoms, Atoms)
        dct = atoms2dict(atoms)
        self.metadata["primitive"] = dct
        talk(".. primitive added to metadata.")

    @property
    def supercell(self):
        """Return the supercell if it is there"""
        if not self._supercell:
            if "supercell" in self.metadata:
                dct = self.metadata["supercell"]
                if "atoms" in dct:
                    dct = dct["atoms"]
                self._supercell = dict2atoms(dct)
            else:
                if self.debug:
                    warn("supercell not provided in trajectory metadata")
        return self._supercell

    @supercell.setter
    def supercell(self, atoms):
        """Set the supercell atoms object"""
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

    @lazy_property
    def cells(self):
        """return cell per time step"""
        return [a.cell[:] for a in self]

    @lazy_property
    def volumes(self):
        """return volume per time step"""
        if all(self.reference_atoms.pbc):
            return [a.get_volume() for a in self]
        return [np.nan for _ in self]

    @property
    def volume(self):
        """return averaged volume"""
        return np.mean(self.volumes).squeeze()

    @property
    def times(self):
        """return the times as numpy array in fs"""
        if self._times is None:
            try:
                fs = self.metadata["MD"]["fs"]
            except KeyError:
                warn("time unit not found in trajectory metadata, use ase.units.fs")
                fs = units.fs

            try:
                times = np.array([a.info["nsteps"] * a.info["dt"] / fs for a in self])
            except KeyError:
                warn("no time steps found, return time as index", level=1)
                times = np.arange(len(self))
            self._times = times
        else:
            assert len(self._times) == len(self)

        return self._times

    @times.setter
    def times(self, new_times):
        """set `trajectory.times`"""
        assert np.size(new_times) == len(self)

        self._times = new_times

    @property
    def timestep(self):
        """return the timestep in fs"""
        from vibes.fourier import get_timestep

        return get_timestep(self.times)

    @lazy_property
    def temperatures(self):
        """return the temperatues as 1d array"""
        return np.array([a.get_temperature() for a in self])

    @property
    def aims_uuid(self):
        """return aims uuids as list"""
        if self._uuid is None:
            self._uuid = [a.info.get(keys.aims_uuid) for a in self]
        return self._uuid

    @aims_uuid.setter
    def aims_uuid(self, uuids: list):
        """set aims uuids (when reading back trajectory.nc)"""
        assert len(uuids) == len(self)
        for uuid, atoms in zip(uuids, self):
            atoms.info[keys.aims_uuid] = str(uuid)

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
        velocities = np.array([a.get_velocities() for a in self])
        if None in velocities:
            velocities = np.full_like(self.positions, np.nan)
        return velocities

    @lazy_property
    def momenta(self):
        """return the velocities as [N_t, N_a, 3] array"""
        return np.array([a.get_momenta() for a in self])

    @lazy_property
    def forces(self):
        """return the forces as [N_t, N_a, 3] array"""
        return np.array([a.get_forces() for a in self])

    @property
    def force_constants_raw(self):
        """return (reduced) force constants or warn if not set"""
        fc = self.metadata[keys.fc]
        if any(x is None for x in (fc, self.primitive, self.supercell)):
            if self.debug:
                warn("`trajectory.force_constants` not set, return None")
        else:
            fc = np.asarray(fc)
            Np, Na = len(self.primitive), len(self.supercell)
            assert fc.shape == (Np, Na, 3, 3), fc.shape

            return fc

    def set_force_constants(self, fc=None):
        """attach force constants as ForceConstants object"""
        if fc is None:
            fc = self.force_constants_raw
        fcs = ForceConstants(
            force_constants=fc, primitive=self.primitive, supercell=self.supercell
        )
        self._force_constants = fcs
        self.set_forces_harmonic()

    @property
    def force_constants(self):
        """return ForceConstants object representing the force constants"""
        return self._force_constants

    @property
    def force_constants_remapped(self):
        """return remapped force constants [3 * Na, 3 * Na]"""
        return self.force_constants.remapped

    def set_forces_harmonic(self):
        """compute harmonic force computed from self.force_constants"""
        timer = Timer("Set harmonic forces")
        displacements = self.displacements
        force_constants = self.force_constants_remapped
        forces_ha = [-force_constants @ d.flatten() for d in displacements]
        self._forces_harmonic = np.array(forces_ha).reshape(self.positions.shape)
        timer()

    @property
    def forces_harmonic(self):
        """return harmonic forces, None if not set via `set_force_constants`"""
        return self._forces_harmonic

    @lazy_property
    def kinetic_energy(self):
        """return the kinetic energy as [N_t] array"""
        return np.array([a.get_kinetic_energy() for a in self])

    @lazy_property
    def potential_energy(self):
        """return the potential energy as [N_t] array"""
        return np.array([a.get_potential_energy() for a in self])

    @property
    def potential_energy_harmonic(self):
        if self._forces_harmonic is not None:
            return -0.5 * (self.forces_harmonic * self.displacements).sum(axis=(1, 2))

    @lazy_property
    def stress_potential(self):
        """return the potential stress as [N_t, 3, 3] array"""
        zeros = np.zeros((3, 3))
        stresses = []
        timer = Timer("Get potential stress")
        for a in progressbar(self):
            if has_stress(a):
                stress = a.get_stress(voigt=False, include_ideal_gas=False)
            else:
                stress = np.full_like(zeros, np.nan)

            stresses.append(stress)
        timer()

        return np.array(stresses, dtype=float)

    @lazy_property
    def stress_kinetic(self):
        """return the kinetic stress as [N_t, 3, 3] array"""
        zeros = np.zeros((3, 3))
        stresses = []
        timer = Timer("Get kinetic stress")
        for ii, a in enumerate(progressbar(self)):
            if has_stress(a):
                stress_potential = self.stress_potential[ii]
                stress_full = a.get_stress(voigt=False, include_ideal_gas=True)
                stress_kinetic = stress_full - stress_potential
            else:
                stress_kinetic = np.full_like(zeros, np.nan)
            stresses.append(stress_kinetic)
        timer()

        return np.array(stresses, dtype=float)

    @property
    def stress(self):
        """return the full stress (kinetic + potential) as [N_t, 3, 3] array"""
        return self.stress_potential + self.stress_kinetic

    @lazy_property
    def stresses_potential(self):
        """return the atomic stress as [N_t, N_a, 3, 3] array"""
        atomic_stresses = []

        zeros = np.zeros((len(self.reference_atoms), 3, 3))

        timer = Timer("Get potential stresses")
        for a in progressbar(self):
            if has_stresses(a):
                atomic_stress = get_stresses(a)
            else:
                atomic_stress = np.full_like(zeros, np.nan)
            atomic_stresses.append(atomic_stress)
        timer()

        if not np.isnan(atomic_stresses).all():
            return np.array(atomic_stresses, dtype=float)

    @lazy_property
    def virials(self):
        """return the virials as [N_t, N_a, 3, 3] array"""
        virials = []

        zeros = np.zeros((len(self.reference_atoms), 3, 3))

        timer = Timer("Get virials")
        for a in self:
            if has_virials(a):
                atom_virials = get_virials(a)
            else:
                atom_virials = np.full_like(zeros, np.nan)
            virials.append(atom_virials)
        timer()

        if not np.isnan(virials).all():
            return np.array(virials, dtype=float)

    def get_pressure(self, kinetic=False):
        """return the pressure as [N_t] array"""
        if kinetic:
            stress = self.stress_kinetic
        else:
            stress = self.stress_potential

        pressure = np.array([-1 / 3 * np.trace(s) for s in stress])
        assert len(pressure) == len(self)

        return pressure  # clean_pressure(pressure)

    @lazy_property
    def pressure_potential(self):
        """return the potential pressure as [N_t] array"""
        return self.get_pressure(kinetic=False)

    @lazy_property
    def pressure_kinetic(self):
        """return the potential pressure as [N_t] array"""
        return self.get_pressure(kinetic=True)

    @lazy_property
    def pressure(self):
        """return the full pressure (kinetic + potential) as [N_t] array"""
        return self.pressure_kinetic + self.pressure_potential

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

    @property
    def dataset(self):
        """return data as xarray.Dataset

        Contains:
            positions, velocities, forces, stress, pressure, temperature
        """
        from .dataset import get_trajectory_dataset

        return get_trajectory_dataset(self)

    @property
    def dataframe(self):
        """return 1D data as pandas.DataFrame

        Contains:
            temperature, kinetic energy, potential energy, pressure
        """
        _keys = [
            keys.temperature,
            keys.energy_kinetic,
            keys.energy_potential,
            keys.pressure_kinetic,
            keys.pressure_potential,
            keys.pressure,
        ]
        df = self.dataset[_keys].to_dataframe()

        return df

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
        """Clean constant drift CAUTION: respect ASE time unit correctly!"""

        timer = Timer("Clean trajectory from constant drift")

        p_drift = np.mean([a.get_momenta().sum(axis=0) for a in self], axis=0)

        talk(f".. drift momentum is {p_drift}")

        for atoms, time in zip(self, self.times):
            atoms.set_momenta(atoms.get_momenta() - p_drift / len(atoms))

            # the displacement
            disp = p_drift / atoms.get_masses().sum() * time
            atoms.positions = atoms.positions - disp

        timer("velocities and positions cleaned from drift")

    def write(self, file=filenames.trajectory):
        """Write to son or nc file

        Args:
            file: path to trajecotry son or nc file
        """
        from .io import write

        write(self, file=file)

    def to_xyz(self, file="positions.xyz", extxyz=False):
        """Write positions to simple xyz file for e.g. viewing with VMD

        Args:
            file: path to trajecotry xyz file
        """
        if extxyz:
            from ase.io.extxyz import write_extxyz as write

            kw = {"write_results": False}
        else:
            from ase.io.xyz import simple_write_xyz as write

            kw = {}

        with open(file, "w") as fo:
            write(fo, self, **kw)

    def to_tdep(self, folder=".", skip=1, stride=1):
        """Convert to TDEP infiles for direct processing

        Args:
            folder: Directory to store tdep files
            skip: Number of structures to skip
        """
        from .io import to_tdep

        to_tdep(self, folder, skip, stride)

    def to_db(self, database):
        """Convert to ase database

        Args:
            database: Filename or address of database

        """
        from .io import to_db

        to_db(self, database)

    def to_traj(self, trajectory):
        """Convert to ase trajectory

        Args:
            trajectory: Filename of ase trajectory

        """
        from .io import to_ase_trajectory

        to_ase_trajectory(self, trajectory)

    def set_displacements(self):
        """calculate the displacements for `reference_atoms`"""
        if not self.supercell:
            # warn("Supercell not set, let us stop here.", level=2)
            warn("SUPERCELL NOT SET, compute w.r.t to reference atoms", level=1)

        timer = Timer("Compute displacements")

        cell = np.asarray(self.reference_atoms.cell)
        shape = self.positions.shape
        displacements = self.positions - self.reference_atoms.positions

        displacements = find_mic(displacements.reshape(-1, 3), cell)[0]
        self.displacements = displacements.reshape(*shape)
        timer()

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
        """Return averaged positions

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

        timer = Timer("Get Atoms hashes")
        hashes = []
        try:
            hashes = [None for aa in range(self[-1].info["displacement_id"] + 1)]
            for atoms in self:
                hashes[atoms.info["displacement_id"]] = hash_atoms(atoms)
        except KeyError:
            for atoms in self:
                hashes.append(hash_atoms(atoms))
        timer()

        return hashes

    # @lazy_property
    # def hash(self):
    #     """hash the atoms and metadata"""
    #     hashes = self.get_hashes()

    #     return hashfunc("".join(hashes))

    @property
    def hash_raw(self):
        """raw hash for input trajectory file (like trajectory.son)"""
        return self._hash_raw

    @hash_raw.setter
    def hash_raw(self, string: str):
        """set raw hash for input trajectory file, can be set only once"""
        assert self._hash_raw is None
        assert isinstance(string, str)
        self._hash_raw = string

    def compute_heat_flux(self):
        """attach `heat_flux` to each `atoms`"""

        have_virials = self.virials is not None
        have_stresses = self.stresses_potential is not None

        if have_virials:
            if have_stresses:
                talk(f"We have both virials and stresses, using virials.")

            talk(f"Computing heat flux using virials.")
            raw_virials = self.virials
        else:
            if not have_stresses:
                talk(f"We have neither virials nor stresses, can't compute heat flux.")
                return None

            talk(f"Computing heat flux using stresses.")
            raw_virials = self.stresses_potential

        # 0) obtain virials and decide for which timesteps to compute heat flux
        idx = []
        all_virials = []
        for i, virials in enumerate(raw_virials):
            if not np.isnan(virials).any():
                idx.append(i)
                all_virials.append(virials)

        # 1) compute average virials
        avg_virials = np.mean(all_virials, axis=0)

        # 2) compute J_avg from average virials
        timer = Timer("Compute heat flux:")
        for i, ii in progressbar(list(enumerate(idx))):
            a = self[ii]
            virials = all_virials[i]
            ds = virials - avg_virials

            # velocity in \AA / fs
            vs = a.get_velocities() * units.fs  # * 1000

            fluxes = np.squeeze(ds @ vs[:, :, None])
            fluxes_aux = np.squeeze(avg_virials @ vs[:, :, None])

            flux = fluxes.sum(axis=0)
            flux_aux = fluxes_aux.sum(axis=0)

            d = {
                keys.heat_flux: flux,
                keys.heat_flux_aux: flux_aux,
            }

            a.calc.results.update(d)

        # attach attributes
        self.set_heat_flux()

        timer()

    def set_heat_flux(self) -> None:
        flux = []
        flux_aux = []
        nan = np.full_like(np.zeros(3), np.nan)
        # attach heat_flux and heat_flux_aux
        for a in self:
            try:
                f = a.calc.results[keys.heat_flux]
            except KeyError:
                f = nan
            flux.append(f)
            try:
                f = a.calc.results[keys.heat_flux_aux]
            except KeyError:
                f = nan
            flux_aux.append(f)
        # set attributes if not all are None
        if not np.isnan(flux).all():
            self._heat_flux = np.array(flux, dtype=float)
            self._heat_flux_aux = np.array(flux_aux, dtype=float)

    def get_heat_flux(self, aux=False):
        """return the heat flux as [N_t, 3] array"""
        # return
        if aux:
            flux = self._heat_flux_aux
        else:
            flux = self._heat_flux

        return flux

    def compute_heat_flux_harmonic(self):
        """compute harmonic properties from FCCalculator"""
        from vibes.ase.calculators.fc import FCCalculator

        calc_ha = FCCalculator(
            atoms_reference=self.reference_atoms,
            force_constants=self.force_constants.remapped,
        )

        timer = Timer("Set harmonic properties")
        for atoms in progressbar(self):
            calc_ha.calculate(atoms)
            for key in calc_ha.results:
                atoms.calc.results[f"{key}_harmonic"] = calc_ha.results[key]

        # set attribute
        self.set_heat_flux_harmonic()

        timer()

    def set_heat_flux_harmonic(self) -> None:
        """attach heat_flux_harmonic as attribute"""
        flux = [a.calc.results[keys.heat_flux_harmonic] for a in self]
        self._heat_flux_harmonic = np.array(flux)
        flux = [a.calc.results[keys.heat_flux_0_harmonic] for a in self]
        self._heat_flux_0_harmonic = np.array(flux)

    def get_heat_flux_harmonic(self) -> np.ndarray:
        """return the harmonic heat flux as [N_t, 3] array"""
        return self._heat_flux_harmonic

    def get_heat_flux_0_harmonic(self) -> np.ndarray:
        """return the harmonic heat flux w.r.t. fixed positions as [N_t, 3] array"""
        return self._heat_flux_0_harmonic
