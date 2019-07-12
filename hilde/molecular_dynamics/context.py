"""Phonopy workflow context managing"""

from pathlib import Path

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

    @property
    def maxsteps(self):
        """return the maxsteps from settings"""
        return self.settings.obj["maxsteps"]
