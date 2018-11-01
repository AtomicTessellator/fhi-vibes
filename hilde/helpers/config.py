import time
import configparser
import json
from collections import OrderedDict

class ClassDict(OrderedDict):
    """Dictionary that supports .attribute access."""
    def __init__(self, *args, **kwargs):
        super(ClassDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class Config(configparser.ConfigParser):
    """ConfigParser that has a slightly more clever get function."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs,
                         interpolation=configparser.ExtendedInterpolation())

    def getval(self, *args, **kwargs):
        try:
            return json.loads(self.get(*args, **kwargs))
        except json.JSONDecodeError:
            try:
                return self.getboolean(*args, **kwargs)
            except ValueError:
                return self.get(*args, **kwargs)

class ConfigDict(ClassDict):
    """Dictionary that holds the configuration settings"""
    def __init__(self, config_files=['highaims.cfg'],
                 *args, **kwargs):

        super().__init__(*args, **kwargs)

        config = Config()
        config.read(config_files)

        # Recursion depth: 1
        for sec in config.sections():
            self[sec] = ClassDict()
            for key in config[sec]:
                val = config.getval(sec, key)
                self[sec][key] = val

    def write(self, filname='config.ini', pickle=False):
        """write a settings object human readable and pickled"""
        with open(filname, 'w') as f:
            timestr = time.strftime('%Y/%m/%d %H:%M:%S')
            f.write(f'# configfile written at {timestr}\n')
            f.write(self.get_string())
        #
        if pickle:
            import pickle
            # write pickled
            with open(filname + '.pick', 'wb') as f:
                pickle.dump(self, f)

    def get_string(self):
        string = ''
        for sec in self:
            string += f'\n[{sec}]\n'
            for key in self[sec]:
                elem = self[sec][key]
                if 'numpy.ndarray' in str(type(elem)):
                    elem = elem.flatten()
                #
                if key == 'verbose':
                    continue
                string += '{:20s} {}\n'.format(f'{key}:', elem)
        return string
