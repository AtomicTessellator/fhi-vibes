"""Phonopy workflow context managing"""

from pathlib import Path

from ase.io import read
from hilde.settings import TaskSettings
from ._defaults import defaults, name, mandatory_base, mandatory_task


class MDSettings(TaskSettings):
    """MD settings. Ensures that settings.md is set up sensibly"""

    def __init__(self, settings):
        """Settings in the context of an md workflow

        Parameters
        -----------
        settings: Settings
            name of the settings file (phonopy.in)
        """

        super().__init__(
            name,
            settings=settings,
            defaults=defaults,
            mandatory_keys=mandatory_base,
            mandatory_obj_keys=mandatory_task,
        )


class MDContext:
    """context for phonopy calculation"""

    def __init__(self, settings, workdir=None, trajectory=None):
        """Initialization

        Parameters
        ----------
        settings: Settings
            settings for the MD Workflow
        workdir: str or Path
            Working directory for the MD workflow
        trajectory: str or Path
            Path to output trajectory
        """
        self.settings = MDSettings(settings)

        self.workdir = self.settings.workdir

        if workdir:
            self.workdir = Path(workdir)
        if not self.workdir:
            self.workdir = name

        if trajectory:
            self.trajectory = Path(trajectory)
        else:
            self.trajectory = Path(self.workdir) / "trajectory.son"

        self._primitive = None
        self._supercell = None

    @property
    def maxsteps(self):
        """return the maxsteps from settings"""
        return self.settings.obj["maxsteps"]

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
