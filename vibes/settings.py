""" Settings class for holding settings, based on configparser.ConfigParser """
from pathlib import Path

from jconfigparser import Config

from vibes._defaults import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_FIREWORKS_FILE,
    DEFAULT_GEOMETRY_FILE,
    DEFAULT_SETTINGS_FILE,
)
from vibes.helpers.attribute_dict import AttributeDict
from vibes.helpers.warnings import warn


class SettingsError(Exception):
    """error in settings"""


def verify_key(
    key: str,
    obj: dict,
    hint: str = None,
    section: bool = False,
    allowed_to_fail: bool = False,
):
    """verify that key is in object, otherwise raise SettingsError

    Args:
        key: Key to check if it is in obj
        obj: Dict to see if key is in it
        hint: string representation of obj
        section: If True key is a section in obj
        allowed_to_fail: If True use wannings not errors
    """
    if not hint:
        hint = str(obj)

    if key not in obj:
        if section:
            msg = f"\n  section [{key}] is missing in {hint}"
        else:
            msg = f"\n  key '{key}' is missing in {hint}"

        if allowed_to_fail:
            warn(msg, level=1)
        else:
            raise SettingsError(msg)


class Configuration(Config):
    """class to hold the configuration from .vibesrc"""

    def __init__(self, config_file: str = DEFAULT_CONFIG_FILE):
        """Initializer

        Args:
            config_file: Path to the configure file
        """
        from vibes import __version__ as version

        super().__init__(filenames=config_file)

        # include the vibes version tag
        self.update({"vibes": {"version": version}})


class Settings(Config):
    """Class to hold the settings parsed from settings.in (+ the configuration)"""

    def __init__(
        self,
        settings_file: str = DEFAULT_SETTINGS_FILE,
        read_config: bool = True,
        config_file: str = DEFAULT_CONFIG_FILE,
        fireworks_file: str = DEFAULT_FIREWORKS_FILE,
        dct: dict = None,
    ):
        """Initialize Settings

        Args:
            settings_file: Path to the settings file
            read_config: read the configuration files
            config_file: Path to the configuration file
            fireworks_file: Path to the FireWorks Configuration file
            dct: create Settings from this dictionary
        """
        if read_config:
            config_files = [config_file, settings_file, fireworks_file]
        else:
            config_files = [settings_file]

        if dct:
            # recursion depth 1
            super().__init__()
            for key in dct:
                self[key] = dct[key]
            if hasattr(dct, "settings_file"):
                self._settings_file = dct.settings_file
        else:
            super().__init__(filenames=[file for file in config_files if file])

        self._settings_file = settings_file

    @classmethod
    def from_dict(cls, dct):
        """initialize from dictionary"""
        return cls(dct=dct)

    @property
    def file(self):
        """return path to the settings file"""
        return self._settings_file

    @property
    def settings_file(self):
        """return path to the settings file"""
        return self._settings_file

    def write(self, file=None):
        """write settings to file"""

        if not file:
            file = self.settings_file

        if not Path(file).exists():
            super().write(file=file)
        else:
            warn(f"{file} exists, do not overwrite settings.", level=1)

    def print(self, only_settings: bool = False):
        ignore_sections = None
        if only_settings:
            ignore_sections = Configuration().keys()
        super().print(ignore_sections=ignore_sections)


class SettingsSection(AttributeDict):
    """Wrapper for a section of settings.in"""

    def __init__(self, name, settings=None, defaults=None, mandatory_keys=None):
        """Initialize Settings in a specific context

        Parameters
        ----------
        name: str
            name of the section
        settings: Settings
            Settings object
        defaults: dict
            dictionary with default key/value pairs
        mandatory_keys: list
            mandatory keys in the section
        """

        if defaults is None:
            defaults = {}
        if mandatory_keys is None:
            mandatory_keys = []

        super().__init__(settings[name])

        self._name = name
        self._settings_file = settings.settings_file

        # validate mandatory keys
        for key in mandatory_keys:
            self.verify_key(key)

        for key in defaults.keys():
            self[key] = self.get(key, defaults[key])

    @property
    def name(self):
        """the name of the task/context"""
        return self._name

    def verify_key(self, key):
        """verify that key is in self.obj

        Parameters
        ----------
        key: str
            key to verify is in self.obj

        Raises
        ------
        SettingsError
            If key is not in self.obj
        """
        verify_key(key, self, hint=f"{self._settings_file}, section [{self.name}]")


class TaskSettings(Settings):
    """Wrapper for Settings in the context of a workflow"""

    def __init__(
        self,
        name=None,
        settings=None,
        read_config=True,
        default_kwargs=None,
        mandatory_keys=None,
        obj_key=None,
        mandatory_obj_keys=None,
        debug=False,
    ):
        """Initialize Settings in a specific context

        Parameters
        ----------
        name: str
            name of the context or workflow
        settings_file: str or Path
            location of settings file. Otherwise inferred from name
        read_config: boolean
            read the configuration file, otherwise just use settings
        defaults: dict
            dictionary with default key/value pairs
        mandatory_keys: list
            mandatory keys in `settings`
        mandatory_obj_keys: list
            mandatory keys in `settings.name`

        Attributes
        ----------
        _obj: dict
            this holds the sub dict with name `name`

        """
        if default_kwargs is None:
            default_kwargs = {}
        if mandatory_keys is None:
            mandatory_keys = []
        if mandatory_obj_keys is None:
            mandatory_obj_keys = []
        if settings is None:
            settings = Settings(read_config=read_config)

        # read the bare settings
        super().__init__(dct=settings)

        self._atoms = None
        self._workdir = None
        self._debug = debug
        self._obj = {}

        for key, val in settings.items():
            self[key] = val

        self._name = name

        # verify name
        self.verify_key(name)

        if not obj_key:
            obj_key = name

        # validate mandatory keys
        for key in mandatory_keys:
            self.verify_key(key)

        if obj_key:
            s = SettingsSection(obj_key, settings, default_kwargs, mandatory_obj_keys)
            self[obj_key] = s
            self._obj = self[obj_key]

        # workdir
        if "workdir" in self.obj:
            self.workdir = self.obj.pop("workdir")

        # make sure atoms are read once
        _ = self.atoms

    @property
    def name(self):
        """the name of the task/context"""
        return self._name

    def verify_key(self, key):
        """verify that key is in self

        Parameters
        ----------
        key: str
            section key to check is in self

        Raises
        ------
        SettingsError
            If key is not in self.obj
        """
        verify_key(
            key, self, hint=f"{self.settings_file}", section=True, allowed_to_fail=True
        )

    @property
    def atoms(self):
        """ Return the settings.atoms object """
        if not self._atoms:
            self._atoms = self.get_atoms()
        return self._atoms

    @atoms.setter
    def atoms(self, obj):
        """Set the settings._atoms with an  ase.atoms.Atoms object"""
        from ase.atoms import Atoms

        assert isinstance(obj, Atoms), type(obj)
        self._atoms = obj

    def get_atoms(self, format="aims"):
        """parse the geometry described in settings.in and return as atoms

        Parameters
        ----------
        format: str
            format of self.geometry.file
        """
        from ase.io import read

        # use the file specified in geometry.file or the default (geometry.in)
        if "geometry" in self and "file" in self.geometry and self.geometry.file:
            path = Path(self.geometry.file)
            file = next(path.parent.glob(path.name))
        else:
            file = DEFAULT_GEOMETRY_FILE

        if Path(file).exists():
            return read(file, format=format)

        if self._debug:
            warn(f"Geometry file {file} not found.", level=1)

        return None

    @property
    def obj(self):
        """the object holding the specific settings for the task"""
        return self._obj

    @property
    def workdir(self):
        """wrapper for the working directory"""
        return self._workdir

    @workdir.setter
    def workdir(self, workdir):
        """wrapper for the working directory"""
        self._workdir = workdir
