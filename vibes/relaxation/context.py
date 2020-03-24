"""Phonopy workflow context managing"""

from pathlib import Path

from ase import Atoms, optimize
from ase.calculators.calculator import Calculator
from ase.constraints import ExpCellFilter

from vibes import keys, son
from vibes.aims.context import AimsContext
from vibes.filenames import filenames
from vibes.helpers import talk, warn
from vibes.helpers.converters import input2dict
from vibes.settings import TaskSettings

from . import _defaults as defaults
from .workflow import _prefix, run_relaxation

# duck type ExpCellFilter, see https://gitlab.com/ase/ase/-/issues/603
# is resolved:


class MyExpCellFilter(ExpCellFilter):
    def get_forces(self, apply_constraint=True):
        """overwrite `apply_constraint` from Expcellfilter"""
        return super().get_forces(apply_constraint=apply_constraint)


class RelaxationSettings(TaskSettings):
    """Relaxation settings. Ensures that settings.relaxation is set up sensibly"""

    def __init__(self, settings):
        """Settings in the context of an md workflow"""

        super().__init__(
            defaults.name,
            settings=settings,
            default_kwargs=defaults.kwargs,
            mandatory_keys=defaults.mandatory_base,
            mandatory_obj_keys=defaults.mandatory_task,
        )


class RelaxationContext:
    """context for relaxation"""

    def __init__(self, settings, workdir=None, trajectory_file=None):
        """Initialization

        Args:
            settings: settings for the relaxation Workflow
            workdir: Working directory for the relaxation workflow
            trajectory_file: Path to output trajectory
        """
        self.settings = RelaxationSettings(settings)

        if workdir:
            self.workdir = Path(workdir).absolute()
        if not self.workdir:
            self.workdir = Path(defaults.name).absolute()

        # workdir has to exist
        Path(self.workdir).mkdir(exist_ok=True, parents=True)

        if trajectory_file:
            self.trajectory_file = Path(trajectory_file)
        else:
            self.trajectory_file = Path(self.workdir) / filenames.trajectory

        self._atoms = None
        self._calculator = None
        self._opt = None

        self.kw = {}
        self.exp_cell_filter_kw = {}

        obj = self.settings.obj
        # save paramters that don't got to the optimizer
        kw = {
            "fmax": obj.pop("fmax"),
            "driver": obj.pop("driver"),
            "unit_cell": obj.pop("unit_cell"),
            "decimals": obj.pop("decimals"),
            "fix_symmetry": obj.pop("fix_symmetry"),
            "symprec": obj.pop("symprec"),
        }
        self.kw = kw

        # save parameters that go to the cellfilter
        kw = {
            "hydrostatic_strain": obj.pop("hydrostatic_strain"),
            "constant_volume": obj.pop("constant_volume"),
            "scalar_pressure": obj.pop("scalar_pressure"),
        }
        self.exp_cell_filter_kw = kw

    # views on self.kw
    @property
    def driver(self):
        return self.kw["driver"]

    @property
    def unit_cell(self):
        return self.kw["unit_cell"]

    @property
    def fmax(self):
        return self.kw["fmax"]

    @property
    def decimals(self):
        return self.kw["decimals"]

    @property
    def fix_symmetry(self):
        return self.kw["fix_symmetry"]

    @property
    def symprec(self):
        return self.kw["symprec"]

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
    def calculator(self):
        """the calculator for running the computation"""
        if not self._calculator:
            # create aims from context and make sure forces are computed
            aims_ctx = AimsContext(settings=self.settings, workdir=self.workdir)
            aims_ctx.settings.obj["compute_forces"] = True

            # atomic stresses
            if self.unit_cell:
                aims_ctx.settings.obj["compute_analytical_stress"] = True

            self._calculator = aims_ctx.get_calculator()

        return self._calculator

    @calculator.setter
    def calculator(self, calculator):
        """set the calculator"""
        assert issubclass(calculator.__class__, Calculator)
        self._calculator = calculator

    @property
    def opt(self):
        """the relaxation algorithm"""
        if not self._opt:
            obj = self.settings.obj
            opt_settings = {"logfile": str(Path(self.workdir) / obj.logfile)}
            obj.update(opt_settings)

            if "bfgs" in self.driver.lower():
                opt = optimize.BFGS(atoms=self.opt_atoms, **obj)
            else:
                warn(f"Relaxation driver {obj.driver} not supported.", level=2)

            talk(f"driver: {self.driver}", prefix=_prefix)
            msg = ["settings:", *[f"  {k}: {v}" for k, v in opt.todict().items()]]
            talk(msg, prefix=_prefix)

            self._opt = opt

        return self._opt

    @property
    def opt_atoms(self):
        """return ExpCellFilter(self.atoms, **kwargs) if `unit_cell == True`"""
        kw = self.exp_cell_filter_kw

        if self.fix_symmetry:
            try:
                from ase.spacegroup.symmetrize import FixSymmetry
            except ModuleNotFoundError:
                msg = "`ase.spacegroup.symmetrize.FixSymmetry` is avaible from ASE 3.20"
                raise RuntimeError(msg)

            constr = FixSymmetry(self.atoms, symprec=self.symprec)
            self.atoms.set_constraint(constr)

        if self.unit_cell:
            return MyExpCellFilter(self.atoms, **kw)
        else:
            return self.atoms

    @property
    def metadata(self):
        """return relaxation metadata as dict"""

        opt_dict = self.opt.todict()

        # save time and mass unit
        opt_dict.update(**self.settings.obj)

        # save kws
        opt_dict.update({keys.relaxation_options: self.kw})
        opt_dict.update({keys.expcellfilter: self.exp_cell_filter_kw})

        # other stuff
        dct = input2dict(self.atoms, calculator=self.calculator)

        return {"relaxation": opt_dict, **dct}

    def resume(self):
        """resume from trajectory"""
        return prepare_from_trajectory(self.atoms, self.opt, self.trajectory_file)

    def run(self, timeout=None):
        """run the context workflow"""
        self.resume()
        run_relaxation(self)


def prepare_from_trajectory(atoms, opt, trajectory_file):
    """Take the last step from trajectory and initialize atoms + md accordingly"""

    trajectory_file = Path(trajectory_file)
    if trajectory_file.exists():
        try:
            last_atoms = son.last_from(trajectory_file)
        except IndexError:
            warn(f"** trajectory lacking the first step, please CHECK!", level=2)
        assert "info" in last_atoms["atoms"]
        opt.nsteps = last_atoms["atoms"]["info"]["nsteps"]

        atoms.set_cell(last_atoms["atoms"]["cell"])
        atoms.set_positions(last_atoms["atoms"]["positions"])
        talk(f"Resume relaxation from  {trajectory_file}", prefix=_prefix)
        return True

    talk(f"** {trajectory_file} does not exist, nothing to prepare", prefix=_prefix)
    return False
