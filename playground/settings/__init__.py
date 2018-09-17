from helpers.config import ClassDict, ConfigDict
from playground.helpers.hash import hashfunc 
import time
import pickle

class Settings(ConfigDict):
    """ Class to hold the settings parsed from highaims.conf (or similar)"""
    def __init__(self, config_files='highaims.conf'):
        super().__init__(config_files=config_files)

    def write(self, folname='setup.log'):
        """write a settings object human readable and pickled"""
        with open(folname, 'w') as f:
            timestr = time.strftime('%Y/%m/%d %H:%M:%S')
            f.write(f'# HIGHaims setup from {timestr}\n')
            f.write(self.get_string())
        #
        # write pickled
        with open(folname + '.pick', 'wb') as f:
            pickle.dump(self, f)

    def get_string(self):
        string = ''
        for sec in self:
            string += f'\n[{sec}]\n'
            for key in self[sec]:
                print(self[sec][key])
                elem = self[sec][key]
                if 'numpy.ndarray' in str(type(elem)):
                    elem = elem.flatten()
                #
                if key == 'verbose':
                    continue
                string += '{:20s} {}\n'.format(f'{key}:', elem)
        return string

    def get_hash(self, short = True):
        return hashfunc(self.get_string())

    def for_handler(self, choose=None):
        handler_setup = {**self.machine, **self.watchdog, **self.database, 'k_grid': self.dft.k_grid}
        if choose:
            handler_setup = {**handler_setup, **self[choose]}
        return ClassDict(handler_setup)

    def print(self):
        print(self.get_string())
