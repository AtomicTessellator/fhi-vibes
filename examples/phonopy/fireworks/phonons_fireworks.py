from pathlib import Path

from hilde.settings import Settings
from hilde.fireworks.workflow_generator import generate_workflow

# These line is only used in the example so an absolute path is not used in settings,
# If you set the absolute path in settings you don't need these 2 lines
settings = Settings()
settings.phonopy.workdir = str(Path(settings.phonopy.workdir).absolute())

generate_workflow(settings)
