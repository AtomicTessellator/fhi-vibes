from pathlib import Path

from ase import Atoms
from ase.calculators.calculator import Calculator

from vibes import keys
from vibes.calculator.context import CalculatorContext
from vibes.filenames import filenames
from vibes.settings import Settings


class TaskContext:
    """context for task"""

    def __init__(
        self, settings, name, template_dict=None, workdir=None, trajectory_file=None
    ):
        """Initialization

        Args:
            settings: settings for the task
            name: name of the task
            template_dict: template settings for task
            workdir: Working directory for running the task
            trajectory_file: Path to output trajectory
        """
        self.settings = Settings(dct=settings, template_dict=template_dict)

        self._atoms = None
        self._calculator = None
        self._name = name

        self.kw = self.settings[self._name]

        # workdir has to exist
        if workdir:
            self.kw[keys.workdir] = Path(workdir).absolute()

        Path(self.workdir).mkdir(exist_ok=True, parents=True)

        if trajectory_file:
            self.trajectory_file = Path(trajectory_file)
        else:
            self.trajectory_file = Path(self.workdir) / filenames.trajectory

    @property
    def name(self):
        return self._name

    @property
    def workdir(self):
        """return the working directory"""
        if self.kw[keys.workdir]:
            return Path(self.kw[keys.workdir]).absolute()

    @workdir.setter
    def workdir(self, folder):
        """set and create the working directory. Use a standard name if dir='auto'"""
        self.kw[keys.workdir] = folder
        self.workdir.mkdir(exist_ok=True, parents=True)

    @property
    def atoms(self):
        """The atoms object for running the computation"""
        if not self._atoms:
            self._atoms = self.settings.read_atoms()
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
            aims_ctx = CalculatorContext(
                settings=self.settings, atoms=self.atoms, workdir=self.workdir
            )

            self._calculator = aims_ctx.get_calculator()

        return self._calculator

    @calculator.setter
    def calculator(self, calculator):
        """set the calculator"""
        assert issubclass(calculator.__class__, Calculator)
        self._calculator = calculator
