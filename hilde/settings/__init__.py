from hilde.helpers.config import ClassDict, ConfigDict
from hilde.helpers.hash import hashfunc

class Settings(ConfigDict):
    """ Class to hold the settings parsed from highaims.conf (or similar)"""
    def __init__(self, config_files='highaims.conf'):
        super().__init__(config_files=config_files)

    def get_hash(self, short = True):
        return hashfunc(self.get_string())

    def for_handler(self, choose=None):
        handler_setup = {**self.machine, **self.watchdog, **self.database, 'k_grid': self.dft.k_grid}
        if choose:
            handler_setup = {**handler_setup, **self[choose]}
        return ClassDict(handler_setup)

    def print(self):
        print(self.get_string())
