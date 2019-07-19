"""Phonopy workflow context managing"""

from pathlib import Path

from ase import Atoms, units as u, md as ase_md
from ase.io import read
from ase.calculators.calculator import Calculator
from ase.md.md import MolecularDynamics
from hilde.molecular_dynamics.workflow import run as run_md
from hilde import son
from hilde.aims.context import AimsContext
from hilde.settings import TaskSettings
from hilde.helpers import warn, talk
from ._defaults import defaults, name, mandatory_base, mandatory_task
from .workflow import run_md, _prefix


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

        if workdir:
            self.workdir = Path(workdir).absolute()
        if not self.workdir:
            self.workdir = Path(name).absolute()

        # workdir has to exist
        Path(self.workdir).mkdir(exist_ok=True)

        if trajectory:
            self.trajectory = Path(trajectory)
        else:
            self.trajectory = Path(self.workdir) / "trajectory.son"

        self._atoms = None
        self._calc = None
        self._md = None
        self._primitive = None
        self._supercell = None

    @property
    def workdir(self):
        """return the working directory"""
        if self.settings.workdir:
            return Path(self.settings.workdir).absolute()

    @workdir.setter
    def workdir(self, folder):
        """set the working directory. Use a standard name if dir='auto'"""
        self.settings.workdir = Path(folder)

    @property
    def maxsteps(self):
        """return the maxsteps from settings"""
        return self.settings.obj["maxsteps"]

    @maxsteps.setter
    def maxsteps(self, steps):
        """set maxsteps"""
        self.settings.obj["maxsteps"] = steps

    @property
    def atoms(self):
        """The atoms object for running the computation"""
        if not self._atoms:
            self._atoms = self.settings.get_atoms()
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        """set the atoms"""
        assert isinstance(atoms, Atoms), atoms
        self._atoms = atoms

    @property
    def calc(self):
        """the calculator for running the computation"""
        if not self._calc:
            # create aims from context and make sure forces are computed
            aims_ctx = AimsContext(settings=self.settings, workdir=self.workdir)
            aims_ctx.settings.obj["compute_forces"] = True

            # atomic stresses
            if self.compute_stresses:
                msg = f"Compute atomic stress every {self.compute_stresses} steps"
                talk(msg, prefix=_prefix)
                aims_ctx.settings.obj["compute_heat_flux"] = True

            self._calc = aims_ctx.get_calculator()

        return self._calc

    @calc.setter
    def calc(self, calc):
        """set the calculator"""
        assert issubclass(calc.__class__, Calculator)
        self._calc = calc

    @property
    def md(self):
        """the MD algorithm"""
        if not self._md:
            obj = self.settings.obj
            md_settings = {
                "atoms": self.atoms,
                "timestep": obj.timestep * u.fs,
                "logfile": Path(self.workdir) / "md.log",  # obj.logfile or "md.log",
            }

            if "verlet" in obj.driver.lower():
                md = ase_md.VelocityVerlet(**md_settings)
            elif "langevin" in obj.driver.lower():
                md = ase_md.Langevin(
                    temperature=obj.temperature * u.kB,
                    friction=obj.friction,
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
    def compute_stresses(self):
        """controls `compute_stresses` in `md` section of settings"""

        compute_stresses = 0
        if "compute_stresses" in self.settings.obj:
            # make sure compute_stresses describes a step length
            compute_stresses = self.settings.obj["compute_stresses"]
            if compute_stresses is True:
                compute_stresses = 1
            elif compute_stresses is False:
                compute_stresses = 0
            else:
                compute_stresses = int(compute_stresses)

        return compute_stresses

    def resume(self):
        """resume from trajectory"""
        prepare_from_trajectory(self.atoms, self.md, self.trajectory)

    def run(self, timeout=None):
        """run the context workflow"""
        self.resume()
        run_md(self, timeout=timeout)


def prepare_from_trajectory(atoms, md, trajectory):
    """Take the last step from trajectory and initialize atoms + md accordingly"""

    trajectory = Path(trajectory)
    if trajectory.exists():
        try:
            last_atoms = son.last_from(trajectory)
        except IndexError:
            warn(f"** trajectory lacking the first step, please CHECK!", level=2)
        assert "info" in last_atoms["atoms"]
        md.nsteps = last_atoms["atoms"]["info"]["nsteps"]

        atoms.set_positions(last_atoms["atoms"]["positions"])
        atoms.set_velocities(last_atoms["atoms"]["velocities"])
        talk(f"Resume from step {md.nsteps} in {trajectory}", prefix=_prefix)
        return True

    talk(f"** {trajectory} does not exist, nothing to prepare", prefix=_prefix)
    return False
