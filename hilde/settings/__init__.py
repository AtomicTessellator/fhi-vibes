""" Settings class for holding settings, based on configparser.ConfigParser """
from hilde import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_FIREWORKS_FILE,
    DEFAULT_SETTINGS_FILE,
    DEFAULT_TEMP_SETTINGS_FILE,
)
from hilde.helpers.config import ConfigDict
from hilde.helpers.hash import hashfunc


class Settings(ConfigDict):
    """ Class to hold the settings parsed from highaims.cfg (or similar)"""

    def __init__(
        self,
        settings_file=DEFAULT_SETTINGS_FILE,
        config_file=DEFAULT_CONFIG_FILE,
        fireworks_file=DEFAULT_FIREWORKS_FILE,
        write=False,
    ):
        super().__init__(config_files=[config_file, settings_file, fireworks_file])

        if write:
            self.write(DEFAULT_TEMP_SETTINGS_FILE)

    def get_hash(self):
        return hashfunc(self.get_string())

    def print(self):
        print(self.get_string())
