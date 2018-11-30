""" general things """

from pathlib import Path

supported_tasks = ["relaxation", "phonopy", "md"]

DEFAULT_CONFIG_FILE = Path(__file__).parents[1] / "config" / "hilde.cfg"
DEFAULT_FIREWORKS_FILE = Path(__file__).parents[1] / "config" / "fireworks.cfg"
DEFAULT_GEOMETRY_FILE = "geometry.in"
DEFAULT_SETTINGS_FILE = "settings.in"
DEFAULT_TEMP_SETTINGS_FILE = "temp_settings.in"

from .settings import Settings