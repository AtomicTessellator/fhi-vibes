""" Settings class for holding settings, based on configparser.ConfigParser """
import time
import configparser
import json

from ase.io import read

from hilde._defaults import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_FIREWORKS_FILE,
    DEFAULT_SETTINGS_FILE,
    DEFAULT_TEMP_SETTINGS_FILE,
)
from hilde import __version__ as version
from hilde.helpers.attribute_dict import AttributeDict


class Config(configparser.ConfigParser):
    """ConfigParser that has a slightly more clever get function."""

    def __init__(self, *args, **kwargs):
        super().__init__(
            *args, **kwargs, interpolation=configparser.ExtendedInterpolation()
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

    def __init__(self, config_files=["hilde.cfg"], *args, **kwargs):

        super().__init__(*args, **kwargs)

        config = Config()
        config.read(config_files)

        # Recursion depth: 1
        for sec in config.sections():
            self[sec] = AttributeDict()
            for key in config[sec]:
                val = config.getval(sec, key)
                self[sec][key] = val

    def print(self):
        print(self.get_string())

    def write(self, filename=DEFAULT_SETTINGS_FILE, pickle=False):
        """write a settings object human readable and pickled"""
        with open(filename, "w") as f:
            timestr = time.strftime("%Y/%m/%d %H:%M:%S")
            f.write(f"# configfile written at {timestr}\n")
            f.write(self.get_string())
        #
        if pickle:
            import pickle

            # write pickled
            with open(filename + ".pick", "wb") as f:
                pickle.dump(self, f)

    def get_string(self, width=30):
        """ return string representation for writing etc. """
        string = ""
        for sec in self:
            string += f"\n[{sec}]\n"
            for key in self[sec]:
                elem = self[sec][key]
                if "numpy.ndarray" in str(type(elem)):
                    elem = elem.flatten()
                #
                if elem is None:
                    elem = "null"
                #
                if key == "verbose":
                    continue
                string += "{:{}s} {}\n".format(f"{key}:", width, elem)
        return string


class Configuration(ConfigDict):
    """ class to hold the configuration from hilde.cfg """

    def __init__(self, config_file=DEFAULT_CONFIG_FILE):
        super().__init__(config_files=[config_file])

        # include the hilde version tag
        self.update({"hilde": {"version": version}})


class Settings(ConfigDict):
    """ Class to hold the settings parsed from settings.in (+ the configuration)"""

    def __init__(
        self,
        settings_file=DEFAULT_SETTINGS_FILE,
        config_file=DEFAULT_CONFIG_FILE,
        fireworks_file=DEFAULT_FIREWORKS_FILE,
        write=False,
    ):
        config_files = [config_file, settings_file, fireworks_file]

        super().__init__(config_files=[file for file in config_files if file])

        if write:
            self.write(DEFAULT_TEMP_SETTINGS_FILE)

    def get_atoms(self, format="aims"):
        """ parse the geometry described in settings.in and return as atoms """
        if "atoms" in self:
            return self.atoms
        elif "geometry" in self:
            if "file" in self.geometry:
                self.atoms = read(self.geometry.file, format=format)
                return self.atoms
