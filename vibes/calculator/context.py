"""Aims workflow context managing"""

from pathlib import Path

import click
from ase import Atoms
from ase.io import read

from vibes.settings import SettingsError, TaskSettings

from . import _defaults as defaults


class AimsSettings(TaskSettings):
    """Aims settings. Ensures that settings are set up sensibly"""

    def __init__(self, settings: dict = None):
        """Settings in the context of a phonopy workflow

        Args:
            settings: Workflow settings for the task

        """

        super().__init__(
            defaults.name,
            settings=settings,
            default_kwargs=defaults.kwargs,
            mandatory_keys=defaults.mandatory_base,
            obj_key=defaults.obj_key,
            mandatory_obj_keys=defaults.mandatory_task,
        )

        # basisset
        if "basisset" in self:
            msg = "`basisset.type` is removed in favour of `basissets.default`. Stop"
            raise RuntimeError(msg)

        if defaults.basisset_key not in self:
            msg = f"basisset not specified in {self.settings_file}"
            raise SettingsError(msg)

        # KS_method
        if "ks_method" in self.control:
            self.control["KS_method"] = self.control.pop("ks_method")


class CalculatorContext:
    """context for aims calculation"""

    def __init__(self, settings: dict, workdir: str = None):
        """Constructor

        Args:
            settings: Settings Object for the Workflow
            workdir: Directory to run the calculation in

        """
        self.settings = AimsSettings(settings)
        self.workdir = self.settings.workdir

        if workdir:
            self.workdir = Path(workdir)
        if not self.workdir:
            self.workdir = Path(defaults.name)

        self._ref_atoms = None
        self._primitive = None
        self._supercell = None
        self._atoms_to_calculate = None

    @property
    def geometry_files(self) -> list:
        """The geometry input files"""
        # find geometries
        files = []
        s = self.settings
        if "file" in s.geometry:
            try:
                try:
                    path = next(Path().glob(s.geometry.file))
                except NotImplementedError:
                    path = Path(s.geometry.file)
                assert path.exists()
            except (AssertionError, StopIteration):
                msg = f"Please inspect [geometry] in {self.settings.settings_file}"
                raise click.FileError(s.geometry["file"], msg)
            files.append(path)
            self.settings.geometry["file"] = str(path)

        if "files" in s.geometry:
            paths = sorted(Path().glob(s.geometry.files))
            for path in paths:
                assert path.exists(), path
            files.extend(paths)

        return files

    @property
    def atoms_to_calculate(self) -> list:
        """The atoms that are supposed to be computed"""
        if not self._atoms_to_calculate:
            files = self.geometry_files
            atoms_list = [read(file, format="aims") for file in files]
            self._atoms_to_calculate = atoms_list
        return self._atoms_to_calculate

    @property
    def primitive(self):
        """The primitive cell structure"""
        g = self.settings.geometry
        if not self._primitive and "primitive" in g:
            self._primitive = read(g["primitive"], format="aims")
        return self._primitive

    @property
    def supercell(self):
        """The supercell structure"""
        g = self.settings.geometry
        if not self._supercell and "supercell" in g:
            self._supercell = read(g["supercell"], format="aims")
        return self._supercell

    @property
    def basisset_location(self):
        """Location of the basis set files"""
        loc = self.settings.machine.basissetloc
        return Path(loc)

    @property
    def ref_atoms(self):
        """The reference structure for the calculation"""
        if not self._ref_atoms:
            self._ref_atoms = read(self.geometry_files[0], format="aims")
        return self._ref_atoms

    @ref_atoms.setter
    def ref_atoms(self, atoms: Atoms):
        """ref_atoms setter

        Args:
            atoms: atoms to set ref_atoms

        """
        assert isinstance(atoms, Atoms)
        self._ref_atoms = atoms

    def get_calculator(self):
        """Get the ASE Calculator based on the context"""
        from .setup import setup_aims

        return setup_aims(self)

    @property
    def calculator(self):
        """return ASE calculator based on this context"""
        return self.get_calculator()

    @property
    def name(self):
        """The name of the calculation"""
        return self.settings.name
