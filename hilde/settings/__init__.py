""" Settings class for holding settings, based on configparser.ConfigParser """
from hilde.helpers.config import ConfigDict
from hilde.helpers.hash import hashfunc


class Settings(ConfigDict):
    """ Class to hold the settings parsed from highaims.cfg (or similar)"""

    def __init__(self, config_files="highaims.cfg"):
        super().__init__(config_files=config_files)

    def get_hash(self):
        return hashfunc(self.get_string())

    def print(self):
        print(self.get_string())
