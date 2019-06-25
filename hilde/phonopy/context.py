"""Phonopy workflow context managing"""

from pathlib import Path

from ase import Atoms
from hilde.settings import WorkflowSettings
from hilde.helpers.numerics import get_3x3_matrix
from ._defaults import defaults, name, mandatory_base, mandatory_task


class PhonopySettings(WorkflowSettings):
    """Phonopy settings. Ensures that settings.phonopy is set up sensibly"""

    def __init__(self, settings):
        """Settings in the context of a phonopy workflow

        Args:
            settings (Settings): Settings object with settings for phonopy
        """
        super().__init__(
            name,
            settings=settings,
            defaults=defaults,
            mandatory_keys=mandatory_base,
            mandatory_obj_keys=mandatory_task,
        )

        # validate
        self.obj["supercell_matrix"] = get_3x3_matrix(self._obj.supercell_matrix)

    @property
    def supercell_matrix(self):
        """return settings.phonopy.supercell_matrix"""
        return self.obj["supercell_matrix"]


class PhonopyContext:
    """context for phonopy calculation"""

    def __init__(self, settings, workdir=None):
        self.settings = PhonopySettings(settings)
        if workdir:
            self.workdir = Path(workdir)
        else:
            if "workdir" in self.settings.obj:
                self.workdir = Path(self.settings.obj.pop("workdir"))
            else:
                self.workdir = Path(name)

        self._ref_atoms = None

    @property
    def q_mesh(self):
        """return the q_mesh from settings"""
        return self.settings.obj["q_mesh"]

    @property
    def ref_atoms(self):
        """return the reference Atoms object for the given context"""
        if not self._ref_atoms:
            self._ref_atoms = self.settings.atoms

        return self._ref_atoms

    @ref_atoms.setter
    def ref_atoms(self, atoms):
        assert isinstance(atoms, Atoms)
        self._ref_atoms = atoms

    @property
    def name(self):
        """return name of the workflow"""
        return self.settings.name
