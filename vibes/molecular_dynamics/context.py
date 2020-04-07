"""Phonopy workflow context managing"""

from pathlib import Path

from ase import md as ase_md
from ase import units as u
from ase.io import read
from ase.md.md import MolecularDynamics

from vibes import son
from vibes.context import TaskContext
from vibes.helpers import talk, warn
from vibes.helpers.converters import input2dict

from ._defaults import name, settings_dict
from .workflow import _prefix, run_md


class MDContext(TaskContext):
    """context for phonopy calculation"""

    def __init__(self, settings=None, workdir=None, trajectory_file=None):
        """Context for MD

        Args:
            settings: settings for the MD Workflow
            workdir: Working directory for the MD workflow
            trajectory_file: Path to output trajectory
        """
        super().__init__(settings, name, template_dict=settings_dict)

        self._md = None
        self._primitive = None
        self._supercell = None

        # legacy
        for key in ("temperature", "friction"):
            if key in self.kw:
                self.kw["kwargs"][key] = self.kw.pop(key)

    @property
    def maxsteps(self):
        """return the maxsteps from settings"""
        return self.kw["maxsteps"]

    @maxsteps.setter
    def maxsteps(self, steps):
        """set maxsteps"""
        self.kw["maxsteps"] = steps

    @property
    def md(self):
        """the MD algorithm"""
        self.mkdir()

        if not self._md:
            obj = self.kw
            md_settings = {
                "atoms": self.atoms,
                "timestep": obj.timestep * u.fs,
                "logfile": Path(self.workdir) / obj["kwargs"].pop("logfile", "md.log"),
            }
            if "verlet" in obj.driver.lower():
                md = ase_md.VelocityVerlet(**md_settings)
            elif "langevin" in obj.driver.lower():
                md = ase_md.Langevin(
                    temperature=obj.kwargs["temperature"] * u.kB,
                    friction=obj.kwargs["friction"],
                    **md_settings,
                )
            else:
                warn(f"MD driver {obj.driver} not supported.", level=2)

            talk(f"driver: {obj.driver}", prefix=_prefix)
            msg = ["settings:", *[f"  {k}: {v}" for k, v in md.todict().items()]]
            talk(msg, prefix=_prefix)

            self._md = md

        return self._md

    @md.setter
    def md(self, md):
        """set the md class"""
        assert issubclass(md.__class__, MolecularDynamics)
        self._md = md

    @property
    def primitive(self):
        """The primitive cell structure"""
        if not self._primitive and "files" in self.settings:
            g = self.settings.files
            if "primitive" in g:
                self._primitive = read(g["primitive"], format="aims")
        return self._primitive

    @property
    def supercell(self):
        """The supercell structure"""
        if not self._supercell and "files" in self.settings:
            g = self.settings.files
            if "supercell" in g:
                self._supercell = read(g["supercell"], format="aims")
        return self._supercell

    @property
    def compute_stresses(self):
        """controls `compute_stresses` in `md` section of settings"""

        compute_stresses = 0
        if "compute_stresses" in self.kw:
            # make sure compute_stresses describes a step length
            compute_stresses = self.kw["compute_stresses"]
            if compute_stresses is True:
                compute_stresses = 1
            elif compute_stresses is False:
                compute_stresses = 0
            else:
                compute_stresses = int(compute_stresses)

        return compute_stresses

    @property
    def metadata(self):
        """return MD metadata as dict"""

        md_dict = self.md.todict()

        # save time and mass unit
        md_dict.update({"fs": u.fs, "kB": u.kB, "dt": self.md.dt, "kg": u.kg})

        # other stuff
        dct = input2dict(
            self.atoms,
            calculator=self.calculator,
            primitive=self.primitive,
            supercell=self.supercell,
        )

        return {"MD": md_dict, **dct}

    def resume(self):
        """resume from trajectory"""
        prepare_from_trajectory(self.atoms, self.md, self.trajectory_file)

    def run(self, timeout=None):
        """run the context workflow"""
        self.resume()
        run_md(self, timeout=timeout)


def prepare_from_trajectory(atoms, md, trajectory_file):
    """Take the last step from trajectory and initialize atoms + md accordingly"""

    trajectory_file = Path(trajectory_file)
    if trajectory_file.exists():
        try:
            last_atoms = son.last_from(trajectory_file)
        except IndexError:
            warn(f"** trajectory lacking the first step, please CHECK!", level=2)
        assert "info" in last_atoms["atoms"]
        md.nsteps = last_atoms["atoms"]["info"]["nsteps"]

        atoms.set_positions(last_atoms["atoms"]["positions"])
        atoms.set_velocities(last_atoms["atoms"]["velocities"])
        talk(f"Resume from step {md.nsteps} in {trajectory_file}", prefix=_prefix)
        return True

    talk(f"** {trajectory_file} does not exist, nothing to prepare", prefix=_prefix)
    return False
