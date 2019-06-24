"""Phonopy workflow context managing"""

from pathlib import Path
import attr

from hilde.settings import WorkflowSettings
from ._defaults import defaults, name, mandatory_base, mandatory_task


class MDSettings(WorkflowSettings):
    """MD settings. Ensures that settings.md is set up sensibly"""

    def __init__(self, settings_file=name + ".in"):
        """Settings in the context of an md workflow

        Args:
            settings_file (str/Path, optional): name of the settings file (phonopy.in)

        """

        super().__init__(
            name,
            settings_file=settings_file,
            defaults=defaults,
            mandatory_keys=mandatory_base,
            mandatory_obj_keys=mandatory_task,
        )


@attr.s
class MDContext:
    """context for phonopy calculation"""

    settings_file = attr.ib(default=name + ".in")
    settings = attr.ib()
    workdir = attr.ib()
    trajectory = attr.ib()

    @settings.default
    def set_settings(self):
        """initialize self.settings"""
        return MDSettings(self.settings_file)

    @workdir.default
    def set_workdir(self):
        """return workdir from settings object"""
        if "workdir" in self.settings.obj:
            return Path(self.settings.obj.pop("workdir"))
        return Path(name)

    @trajectory.default
    def set_trajectory(self):
        """return trajectory path"""
        return Path(self.workdir) / "trajectory.son"

    @property
    def maxsteps(self):
        """return the maxsteps from settings"""
        return self.settings.obj["maxsteps"]

