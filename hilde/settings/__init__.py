from hilde.helpers.config import AttributeDict, ConfigDict
from hilde.helpers.hash import hashfunc

class Settings(ConfigDict):
    """ Class to hold the settings parsed from highaims.cfg (or similar)"""
    def __init__(self, config_files='highaims.cfg'):
        super().__init__(config_files=config_files)
        self.fw_settings.remote_host = self.make_list(self.fw_settings.remote_host)
        self.fw_settings.remote_config_dir = self.make_list(self.fw_settings.remote_config_dir)
        self.fw_settings.tasks2queue = self.make_list(self.fw_settings.tasks2queue)


    def get_hash(self, short = True):
        return hashfunc(self.get_string())

    def for_handler(self, choose=None):
        handler_setup = {**self.machine, **self.watchdog, **self.database, 'k_grid': self.dft.k_grid}
        if choose:
            handler_setup = {**handler_setup, **self[choose]}
        return AttributeDict(handler_setup)

    def print(self):
        print(self.get_string())

    def make_list(val):
        return val.split(",")