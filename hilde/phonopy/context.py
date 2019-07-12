"""Phonopy workflow context managing"""

from pathlib import Path

from ase import Atoms
from hilde.settings import TaskSettings
from hilde.helpers.numerics import get_3x3_matrix
from hilde.structure.misc import get_sysname
from ._defaults import defaults, name, mandatory_base, mandatory_task


class PhonopySettings(TaskSettings):
    """Phonopy settings. Ensures that settings.phonopy is set up sensibly"""

    def __init__(self, settings):
        """Settings in the context of a phonopy workflow

        Parameters
        ----------
        settings: Settings
            Settings object with settings for phonopy
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
        """Intializer

        Parameters
        ----------
        settings: Settings
            Settings for the Workflow
        workdir: str or Path
            The working directory for the workflow
        """
        self.settings = PhonopySettings(settings)
        self._ref_atoms = None

        if workdir:
            self.workdir = workdir
        if not self.workdir:
            self.workdir = name
        else:
            self.workdir = self.workdir

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
