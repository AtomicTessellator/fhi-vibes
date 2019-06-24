"""Phonopy workflow context managing"""

from pathlib import Path
import attr

from ase.io import read
from hilde.settings import WorkflowSettings
from hilde.helpers.numerics import get_3x3_matrix
from ._defaults import defaults, name, mandatory_base, mandatory_task


class PhonopySettings(WorkflowSettings):
    """Phonopy settings. Ensures that settings.phonopy is set up sensibly"""

    def __init__(self, settings_file=name + ".in", input_settings=None):
        """Settings in the context of a phonopy workflow

        Args:
            settings_file (str/Path, optional): name of the settings file (phonopy.in)

        """
        super().__init__(
            name,
            settings_file=settings_file,
            defaults=defaults,
            mandatory_keys=mandatory_base,
            mandatory_obj_keys=mandatory_task,
            input_settings=input_settings,
        )

        # validate
        self.obj["supercell_matrix"] = get_3x3_matrix(self._obj.supercell_matrix)

    @property
    def supercell_matrix(self):
        """return settings.phonopy.supercell_matrix"""
        return self.obj["supercell_matrix"]


@attr.s
class PhonopyContext:
    """context for phonopy calculation"""

    settings_file = attr.ib(default=name + ".in")
    input_settings = attr.ib(default=None)
    settings = attr.ib()
    workdir = attr.ib()
    _ref_atoms = None

    @settings.default
    def set_settings(self):
        """initialize self.settings"""
        return PhonopySettings(self.settings_file, self.input_settings)

    @workdir.default
    def set_workdir(self):
        """return workdir from settings object"""
        if "workdir" in self.settings.obj:
            return Path(self.settings.obj.pop("workdir"))
        return Path(name)

    @property
    def q_mesh(self):
        """return the q_mesh from settings"""
        return self.settings.obj["q_mesh"]

    @property
    def ref_atoms(self):
        if not self._ref_atoms and not self.input_settings.atoms:
            self._ref_atoms = read(self.settings.geometry.file, format="aims")
        elif not self._ref_atoms:
            self._ref_atoms = self.input_settings.atoms
        return self._ref_atoms

    @ref_atoms.setter
    def ref_atoms(self, atoms):
        assert isinstance(atoms, Atoms)
        self._ref_atoms = atoms
