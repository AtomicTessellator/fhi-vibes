"""Phonopy workflow context managing"""

from pathlib import Path

from ase import Atoms, units as u, optimize
from ase.io import read
from ase.calculators.calculator import Calculator
from ase.md.md import MolecularDynamics

from hilde import son
from hilde.aims.context import AimsContext
from hilde.settings import TaskSettings
from hilde.helpers import warn, talk
from hilde.helpers.converters import input2dict
from ._defaults import defaults, name, mandatory_base, mandatory_task
from .workflow import run_relaxation, _prefix


class RelaxationSettings(TaskSettings):
    """Relaxation settings. Ensures that settings.relaxation is set up sensibly"""

    def __init__(self, settings):
        """Settings in the context of an md workflow"""

        super().__init__(
            name,
            settings=settings,
            defaults=defaults,
            mandatory_keys=mandatory_base,
            mandatory_obj_keys=mandatory_task,
        )


class RelaxationContext:
    """context for relaxation"""

    def __init__(self, settings, workdir=None, trajectory=None):
        """Initialization

        Args:
            settings: settings for the relaxation Workflow
            workdir: Working directory for the relaxation workflow
            trajectory: Path to output trajectory
        """
        self.settings = RelaxationSettings(settings)

        if workdir:
            self.workdir = Path(workdir).absolute()
        if not self.workdir:
            self.workdir = Path(name).absolute()

        # workdir has to exist
        Path(self.workdir).mkdir(exist_ok=True, parents=True)

        if trajectory:
            self.trajectory = Path(trajectory)
        else:
            self.trajectory = Path(self.workdir) / "trajectory.son"

        self._atoms = None
        self._calc = None
        self._opt = None

    @property
    def workdir(self):
        """return the working directory"""
        if self.settings.workdir:
            return Path(self.settings.workdir).absolute()

    @workdir.setter
    def workdir(self, folder):
        """set and create the working directory. Use a standard name if dir='auto'"""
        self.settings.workdir = Path(folder)
        self.workdir.mkdir(exist_ok=True, parents=True)

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
            if self.unit_cell:
                aims_ctx.settings.obj["compute_analytical_stress"] = True

            self._calc = aims_ctx.get_calculator()

        return self._calc

    @calc.setter
    def calc(self, calc):
        """set the calculator"""
        assert issubclass(calc.__class__, Calculator)
        self._calc = calc

    @property
    def opt(self):
        """the relaxation algorithm"""
        if not self._opt:
            obj = self.settings.obj
            opt_settings = {"logfile": str(Path(self.workdir) / obj.logfile)}
            obj.update(opt_settings)

            self.fmax = obj.pop("fmax")
            self.driver = obj.pop("driver")
            self.unit_cell = obj.pop("unit_cell")
            self.decimals = obj.pop("decimals")

            if "bfgs" in self.driver.lower():
                opt = optimize.BFGS(atoms=self.atoms, **obj)
            else:
                warn(f"Relaxation driver {obj.driver} not supported.", level=2)

            talk(f"driver: {self.driver}", prefix=_prefix)
            msg = ["settings:", *[f"  {k}: {v}" for k, v in opt.todict().items()]]
            talk(msg, prefix=_prefix)

            self._opt = opt

        return self._opt

    @property
    def metadata(self):
        """return relaxation metadata as dict"""

        opt_dict = self.opt.todict()

        # save time and mass unit
        opt_dict.update(**self.settings.obj)

        # other stuff
        dct = input2dict(self.atoms, calc=self.calc)

        return {"relaxation": opt_dict, **dct}

    def resume(self):
        """resume from trajectory"""
        return prepare_from_trajectory(self.atoms, self.opt, self.trajectory)

    def run(self, timeout=None):
        """run the context workflow"""
        self.resume()
        run_relaxation(self)


def prepare_from_trajectory(atoms, opt, trajectory):
    """Take the last step from trajectory and initialize atoms + md accordingly"""

    trajectory = Path(trajectory)
    if trajectory.exists():
        try:
            last_atoms = son.last_from(trajectory)
        except IndexError:
            warn(f"** trajectory lacking the first step, please CHECK!", level=2)
        assert "info" in last_atoms["atoms"]
        opt.nsteps = last_atoms["atoms"]["info"]["nsteps"]

        atoms.set_cell(last_atoms["atoms"]["cell"])
        atoms.set_positions(last_atoms["atoms"]["positions"])
        talk(f"Resume relaxation from  {trajectory}", prefix=_prefix)
        return True

    talk(f"** {trajectory} does not exist, nothing to prepare", prefix=_prefix)
    return False
