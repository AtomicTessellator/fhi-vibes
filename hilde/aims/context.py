"""Aims workflow context managing"""

from pathlib import Path
import attr
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

    def __init__(self, settings_file=name + ".in", input_settings=None):
        """Settings in the context of a phonopy workflow

        Args:
            settings_file (str/Path, optional): name of the settings file (aims.in)

        """

        super().__init__(
            name,
            settings_file=settings_file,
            defaults=defaults,
            mandatory_keys=mandatory_base,
            obj_key=obj_key,
            mandatory_obj_keys=mandatory_task,
            input_settings=input_settings,
        )

        # basisset
        if not any(key in self for key in mandatory_basisset):
            msg = f"basisset not specified in {self.settings_file}"
            raise SettingsError(msg)


@attr.s
class AimsContext:
    """context for aims calculation"""

    settings_file = attr.ib(default=name + ".in")
    input_settings = attr.ib(default=None)
    settings = attr.ib()
    workdir = attr.ib()
    _ref_atoms = None
    _primitive = None
    _supercell = None
    _atoms_to_calculate = None

    @settings.default
    def set_settings(self):
        """initialize self.settings"""
        return AimsSettings(self.settings_file, self.input_settings)

    @workdir.default
    def set_workdir(self):
        """return workdir from settings object"""
        return Path(name)

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
                msg = f"Please inspect [geometry] in {self.settings_file}"
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
