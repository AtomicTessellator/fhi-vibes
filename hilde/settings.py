""" Settings class for holding settings, based on configparser.ConfigParser """
import time
import configparser
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.atoms import Atoms


from hilde._defaults import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_FIREWORKS_FILE,
    DEFAULT_SETTINGS_FILE,
    DEFAULT_GEOMETRY_FILE,
)
from hilde import __version__ as version
from hilde.helpers.attribute_dict import AttributeDict, MultiOrderedDict
from hilde.helpers.warnings import warn


class SettingsError(Exception):
    """error in settings"""


def verify_key(key, obj, hint=None, section=False, allowed_to_fail=False):
    """verify that key is in object, otherwise raise SettingsError

    Parameters
    ----------
    key: str
        Key to check if it is in obj
    obj: dict like object
        Dict to see if key is in it
    hint: str
        string representation of obj
    section: bool
        If True key is a section in obj
    allowed_to_fail: bool
        If True use wannings not errors
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


class Config(configparser.ConfigParser):
    """ConfigParser that uses JSON to parse the values instead returning stings"""

    def __init__(self, *args, **kwargs):
        super().__init__(
            *args,
            **kwargs,
            interpolation=configparser.ExtendedInterpolation(),
            strict=False,
        )

    def getval(self, *args, **kwargs):
        """ Redifine getval() to allow for json formated values (not only string) """
        try:
            return json.loads(self.get(*args, **kwargs))
        except json.JSONDecodeError:
            try:
                return self.getboolean(*args, **kwargs)
            except ValueError:
                return self.get(*args, **kwargs)


class ConfigDict(AttributeDict):
    """Dictionary that holds the configuration settings"""

    def __init__(self, config_files=None, dct=None, **kwargs):
        """Initialize ConfigDict

        Args:
            config_files (list of str): A list of configure files to read in
            dct (dict): a dictionary
        """
        super().__init__(**kwargs)

        # initialize from config files
        if config_files:
            config = Config()
            config.read(config_files)

            # Recursion depth: 1
            for sec in config.sections():
                self[sec] = AttributeDict()
                for key in config[sec]:
                    self[sec][key] = config.getval(sec, key)

            # check for `output` to resolve multiple options
            if "control" in config and "output" in config["control"].keys():
                kw = {"dict_type": MultiOrderedDict, "strict": False}
                output_cfg = configparser.ConfigParser(**kw)
                # discard config files to avoid double sections
                files = [file for file in config_files if str(file).endswith(".in")]
                output_cfg.read(files)
                values = output_cfg["control"]["output"].split("\n")
                self["control"]["output"] = values

        # initialize from dictionary
        elif dct:
            # Recursion depth: 1
            for sec in dct:
                self[sec] = AttributeDict()
                for key in dct[sec]:
                    self[sec][key] = dct[sec][key]

        else:
            raise ValueError("Either provide `config_files` or a dict in `dct`")

    def __str__(self):
        """ for printing the object """
        return self.get_string()

    def print(self, only_settings=False):
        """ literally print(self) """
        print(self.get_string(only_settings=only_settings), flush=True)

    def write(self, filename=DEFAULT_SETTINGS_FILE):
        """write a settings object human readable

        Parameters
        ----------
        filename: str
            path use to write the file
        """
        with open(filename, "w") as f:
            timestr = time.strftime("%Y/%m/%d %H:%M:%S")
            f.write(f"# configfile written at {timestr}\n")
            f.write(self.get_string())

    def get_string(self, width=30, only_settings=False):
        """ return string representation for writing etc.

        Parameters
        ----------
        width: int
            The width of the string column to print
        only settings: bool
            If True only print the settings

        Returns
        -------
        string: str
            The string representation of the ConfigDict
        """
        if only_settings:
            ref_dict = Configuration()
        else:
            ref_dict = {}

        string = ""
        for sec in self:
            # Filter out the private attributes
            if sec.startswith("_") or sec in ref_dict:
                continue

            string += f"\n[{sec}]\n"
            for key in self[sec]:
                elem = self[sec][key]
                if "numpy.ndarray" in str(type(elem)):
                    elem = np.array2string(elem.flatten(), separator=",")
                #
                if elem is None:
                    elem = "null"
                #
                if key == "verbose":
                    continue
                # write out `output` keys one by one for readability
                if key == "output":
                    for elem in self[sec][key]:
                        string += "{:{}s} {}\n".format(f"{key}:", width, elem)
                else:
                    string += "{:{}s} {}\n".format(f"{key}:", width, elem)
        return string


class Configuration(ConfigDict):
    """class to hold the configuration from .hilderc"""

    def __init__(self, config_file=DEFAULT_CONFIG_FILE):
        """Initializer

        Parameters
        ----------
        config_file: str
            Path to the configure file
        """
        super().__init__(config_files=config_file)

        # include the hilde version tag
        self.update({"hilde": {"version": version}})


class Settings(ConfigDict):
    """Class to hold the settings parsed from settings.in (+ the configuration)"""

    def __init__(
        self,
        settings_file=DEFAULT_SETTINGS_FILE,
        read_config=True,
        config_file=DEFAULT_CONFIG_FILE,
        fireworks_file=DEFAULT_FIREWORKS_FILE,
        dct=None,
    ):
        """Initialize Settings

        Args:
            settings_file (str): Path to the settings file
            read_config (boolean): read the configuration files
            config_file (str): Path to the configuration file
            fireworks_file (str): Path to the FireWorks Configuration file
            dct (dict): create Settings from this dictionary
        """
        self._settings_file = settings_file

        if read_config:
            config_files = [config_file, settings_file, fireworks_file]
        else:
            config_files = settings_file

        if dct:
            super().__init__(dct=dct)
        else:
            super().__init__(config_files=[file for file in config_files if file])

    @classmethod
    def from_dict(cls, dct):
        """initialize from dictionary"""
        return cls(dct=dct)

    @property
    def settings_file(self):
        """return path to the settings file"""
        return self._settings_file

    def write(self, filename=None):
        """write settings to file"""

        if not filename:
            filename = self.settings_file

        if not Path(filename).exists():
            super().write(filename=filename)
        else:
            warn(f"{filename} exists, do not overwrite settings.", level=1)


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
        defaults=None,
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
        if defaults is None:
            defaults = {}
        if mandatory_keys is None:
            mandatory_keys = []
        if mandatory_obj_keys is None:
            mandatory_obj_keys = []

        # read the bare settings
        super().__init__(settings_file=settings.settings_file, read_config=read_config)

        self._atoms = None
        self._workdir = None
        self._debug = debug
        self._obj = {}

        for key, val in settings.items():
            self[key] = val

        self._name = name

        if not obj_key:
            obj_key = name

        # validate mandatory keys
        for key in mandatory_keys:
            self.verify_key(key)

        if obj_key:
            s = SettingsSection(obj_key, settings, defaults, mandatory_obj_keys)
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
        assert isinstance(obj, Atoms), type(obj)
        self._atoms = obj

    def get_atoms(self, format="aims"):
        """parse the geometry described in settings.in and return as atoms

        Parameters
        ----------
        format: str
            format of self.geometry.file
        """

        # use the file specified in geometry.file or the default (geometry.in)
        if "geometry" in self and "file" in self.geometry and self.geometry.file:
            file = next(Path().glob(self.geometry.file))
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
