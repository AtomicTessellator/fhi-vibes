""" useful things to import """
# flake8: noqa

import pkg_resources

from ._defaults import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_FIREWORKS_FILE,
    DEFAULT_GEOMETRY_FILE,
    DEFAULT_SETTINGS_FILE,
    supported_tasks,
)
from .settings import Configuration, Settings

__version__ = str(pkg_resources.require("fhi-vibes")[0].version)


# from .templates.aims import setup_aims
# from .helpers.restarts import restart
# from .tasks import run
# from .trajectory.io import reader
