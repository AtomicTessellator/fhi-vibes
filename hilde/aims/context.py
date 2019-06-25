"""Aims workflow context managing"""

from pathlib import Path
import click

from ase import Atoms
from ase.io import read

from hilde.settings import WorkflowSettings, SettingsError
from .setup import setup_aims
from ._defaults import (
    name,
    obj_key,
    defaults,
    mandatory_base,
    mandatory_task,
    mandatory_basisset,
)


class AimsSettings(WorkflowSettings):
    """Aims settings. Ensures that settings are set up sensibly"""

    def __init__(self, settings=None):
        """Settings in the context of a phonopy workflow

        Args:
            settings (Settings): Settings object

        """

        super().__init__(
            name,
            settings=settings,
            defaults=defaults,
            mandatory_keys=mandatory_base,
            obj_key=obj_key,
            mandatory_obj_keys=mandatory_task,
        )

        # basisset
        if not any(key in self for key in mandatory_basisset):
            msg = f"basisset not specified in {self.settings_file}"
            raise SettingsError(msg)


class AimsContext:
    """context for aims calculation"""

    def __init__(self, settings, workdir=None):
        self.settings = AimsSettings(settings)

        if workdir:
            self.workdir = Path(workdir)
        else:
            self.workdir = Path(name)

        self._ref_atoms = None
        self._primitive = None
        self._supercell = None
        self._atoms_to_calculate = None

    @property
    def geometry_files(self):
        """return the geometry input"""
        # find geometries
        filenames = []
        s = self.settings
        if "file" in s.geometry:
            try:
                path = next(Path().glob(s.geometry.file))
                assert path.exists()
            except (AssertionError, StopIteration):
                msg = f"Please inspect [geometry] in {self.settings.settings_file}"
                raise click.FileError(s.geometry["file"], msg)
            filenames.append(path)

        if "files" in s.geometry:
            paths = sorted(Path().glob(s.geometry.files))
            for path in paths:
                assert path.exists(), path
            filenames.extend(paths)

        return filenames

    @property
    def atoms_to_calculate(self):
        """return atoms that are supposed to be computed"""
        if not self._atoms_to_calculate:
            filenames = self.geometry_files
            atoms_list = [read(file, format="aims") for file in filenames]
            self._atoms_to_calculate = atoms_list
        return self._atoms_to_calculate

    @property
    def primitive(self):
        """return primitive structure"""
        g = self.settings.geometry
        if not self._primitive and "primitive" in g:
            self._primitive = read(g["primitive"], format="aims")
        return self._primitive

    @property
    def supercell(self):
        """return supercell"""
        g = self.settings.geometry
        if not self._supercell and "supercell" in g:
            self._supercell = read(g["supercell"], format="aims")
        return self._supercell

    @property
    def basisset_location(self):
        loc = self.settings.machine.basissetloc
        return Path(loc)

    @property
    def ref_atoms(self):
        if not self._ref_atoms:
            self._ref_atoms = read(self.geometry_files[0], format="aims")
        return self._ref_atoms

    @ref_atoms.setter
    def ref_atoms(self, atoms):
        assert isinstance(atoms, Atoms)
        self._ref_atoms = atoms

    def get_calculator(self):
        """return and ase aims calculator object based on the context"""
        return setup_aims(self)

    @property
    def name(self):
        return self.settings.name
