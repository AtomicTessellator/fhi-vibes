"""Phonopy workflow context managing"""

import sys
from pathlib import Path

from ase import Atoms

from vibes.helpers.numerics import get_3x3_matrix
from vibes.settings import TaskSettings
from vibes.structure.misc import get_sysname

from . import metadata2dict
from . import postprocess as postprocess
from . import wrapper as backend
from ._defaults import defaults, mandatory, name


class PhonopySettings(TaskSettings):
    """Phonopy settings. Ensures that settings.phonopy is set up sensibly"""

    def __init__(self, settings, name=name, defaults_kw=None, mandatory_kw=None):
        """Settings in the context of a phonopy workflow

        Parameters
        ----------
        settings: Settings
            Settings object with settings for phonopy
        """
        if defaults_kw is None:
            defaults_kw = defaults
        if mandatory_kw is None:
            mandatory_kw = mandatory

        super().__init__(
            name=name, settings=settings, defaults=defaults_kw, **mandatory_kw
        )

        # validate
        self.obj["supercell_matrix"] = get_3x3_matrix(self._obj.supercell_matrix)

    @property
    def supercell_matrix(self):
        """return settings.phonopy.supercell_matrix"""
        return self.obj["supercell_matrix"]


class PhonopyContext:
    """context for phonopy calculation"""

    def __init__(
        self,
        settings,
        workdir=None,
        name=name,
        defaults_kw: dict = None,
        mandatory_kw: dict = None,
    ):
        """Intializer

        Parameters
        ----------
        settings: Settings
            Settings for the Workflow
        workdir: str or Path
            The working directory for the workflow
        """
        self.settings = PhonopySettings(
            settings, name=name, defaults_kw=defaults_kw, mandatory_kw=mandatory_kw
        )
        self._ref_atoms = None

        if workdir:
            self.workdir = workdir
        if not self.workdir:
            self.workdir = name
        else:
            self.workdir = self.workdir

        self.backend = backend
        self.postprocess = postprocess

    @property
    def workdir(self):
        """return the working directory"""
        return self.settings.workdir

    @workdir.setter
    def workdir(self, folder):
        """set the working directory. Use a standard name if dir='auto'"""
        if "auto" in str(folder).lower():
            smatrix = self.settings.obj.supercell_matrix.flatten()
            vol = self.ref_atoms.get_volume()
            sysname = get_sysname(self.ref_atoms)
            rep = "_{}_{}{}{}_{}{}{}_{}{}{}_{:.3f}".format(sysname, *smatrix, vol)
            dirname = name + rep
            self.settings.workdir = Path(dirname)
        else:
            self.settings.workdir = Path(folder)

    @property
    def q_mesh(self):
        """return the q_mesh from settings"""
        return self.settings.obj["q_mesh"]

    @property
    def ref_atoms(self):
        """return the reference Atoms object for the given context"""
        if not self._ref_atoms:
            self._ref_atoms = self.settings.atoms.copy()

        return self._ref_atoms

    @ref_atoms.setter
    def ref_atoms(self, atoms):
        """The setter for ref_atoms, makes sure it's an atoms object indeed"""
        assert isinstance(atoms, Atoms)
        self._ref_atoms = atoms

    @property
    def name(self):
        """return name of the workflow"""
        return self.settings.name

    def preprocess(self):
        return self.backend.preprocess(atoms=self.ref_atoms, **self.settings.obj)

    def bootstrap(self, dry=False):
        """load settings, prepare atoms, calculator, and phonopy object"""

        # Phonopy preprocess
        phonon, supercell, scs = self.preprocess()

        # if calculator not given, create an aims context for this calculation
        if self.settings.atoms and self.settings.atoms.calc:
            calc = self.settings.atoms.calc
        else:
            from vibes import aims

            aims_ctx = aims.AimsContext(settings=self.settings, workdir=self.workdir)
            # set reference structure for aims calculation
            # and make sure forces are computed
            aims_ctx.ref_atoms = supercell
            aims_ctx.settings.obj["compute_forces"] = True

            calc = aims.setup.setup_aims(aims_ctx)

        # save metadata
        metadata = metadata2dict(phonon, calc)

        return {
            "atoms_to_calculate": scs,
            "calculator": calc,
            "metadata": metadata,
            "workdir": self.workdir,
            "settings": self.settings,
            "save_input": True,
            "backup_after_calculation": False,
            "dry": dry,
            **self.settings.obj,
        }

    def run(self, dry=False):
        """run phonopy workflow """
        from vibes.helpers import talk
        from vibes.helpers.restarts import restart
        from vibes.tasks import calculate_socket

        args = self.bootstrap(dry=dry)

        talk(f"Run phonopy workflow in working directory\n  {self.workdir}")

        try:
            self.postprocess(workdir=self.workdir)
            msg = "** Postprocess could be performed from previous calculations. Check"
            msg += f"\n**  {self.workdir}"
            sys.exit(msg)
        except (FileNotFoundError, RuntimeError):
            completed = calculate_socket(**args)

        if not completed:
            restart(args["settings"])
        else:
            talk("Start postprocess.")
            self.postprocess(**args)
            talk("done.")
