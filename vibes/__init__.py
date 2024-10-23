"""useful things to import"""
# flake8: noqa

from importlib.metadata import version 

from ._defaults import (
    DEFAULT_CONFIG_FILE,
    DEFAULT_FIREWORKS_FILE,
    DEFAULT_GEOMETRY_FILE,
    DEFAULT_SETTINGS_FILE,
    supported_tasks,
)
from .settings import Configuration, Settings


__version__ = version("fhi-vibes")
