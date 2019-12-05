from ase.io import read

from fireworks import Workflow
from fireworks.fw_config import LAUNCHPAD_LOC

from pathlib import Path

from vibes import DEFAULT_CONFIG_FILE
from vibes.fireworks.launchpad import LaunchPadHilde
from vibes.settings import Settings
from vibes.fireworks.workflows.workflow_generator import generate_workflow

settings = Settings()
settings.phono3py.workdir = str(Path(settings.phono3py.workdir).absolute())
atoms = read(settings.geometry.file)

launchpad = LaunchPadHilde.from_file(LAUNCHPAD_LOC)
launchpad.add_wf(generate_workflow(settings, atoms))
