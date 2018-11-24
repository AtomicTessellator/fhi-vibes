""" general things """

from pathlib import Path

supported_tasks = ["relaxation", "phonopy", "md"]

DEFAULT_SETTINGS_FILE = "settings.in"
DEFAULT_GEOMETRY_FILE = "geometry.in"
DEFAULT_TEMP_SETTINGS_FILE = "temp_settings.in"
DEFAULT_CONFIG_FILE = Path(__file__).parents[1] / "config" / "hilde.cfg"
