""" Settings class for holding settings, based on configparser.ConfigParser """
from hilde import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_SETTINGS_FILE,
    DEFAULT_TEMP_SETTINGS_FILE,
    DEFAULT_GEOMETRY_FILE,
)
from hilde.helpers.config import ConfigDict
from hilde.helpers.hash import hashfunc


class Configuration(ConfigDict):
    """ class to hold the configuration from hilde.cfg """

    def __init__(self, config_file=DEFAULT_CONFIG_FILE):
        super().__init__(config_files=[config_file])


class Settings(ConfigDict):
    """ Class to hold the settings parsed from settings.in (+ the configuration)"""

    def __init__(
        self,
        settings_file=DEFAULT_SETTINGS_FILE,
        config_file=DEFAULT_CONFIG_FILE,
        write=False,
    ):

        config_files = [config_file, settings_file]

        super().__init__(config_files=[file for file in config_files if file])

        if write:
            self.write(DEFAULT_TEMP_SETTINGS_FILE)
