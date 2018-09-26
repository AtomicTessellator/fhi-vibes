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
    def __init__(self, config_files=['highaims.conf'],
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
