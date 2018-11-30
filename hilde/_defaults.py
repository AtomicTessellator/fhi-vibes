""" some default naming """

from pathlib import Path

supported_tasks = ["relaxation", "phonopy", "md"]

DEFAULT_CONFIG_FILE = Path(__file__).parents[1] / "hilde.cfg"
DEFAULT_FIREWORKS_FILE = Path(__file__).parents[1] / "fireworks.cfg"
DEFAULT_GEOMETRY_FILE = "geometry.in"
DEFAULT_SETTINGS_FILE = "settings.in"
DEFAULT_TEMP_SETTINGS_FILE = "temp_settings.in"
