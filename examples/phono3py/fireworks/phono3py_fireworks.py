from ase.io import read

from fireworks import Workflow

from hilde import DEFAULT_CONFIG_FILE
from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.settings import Settings
from hilde.workflows.workflow_generator import generate_workflow

settings = Settings()
atoms = read(settings.geometry.file)

launchpad = LaunchPadHilde.from_file("/home/purcell/.fireworks/my_launchpad.yaml")
launchpad.add_wf(generate_workflow(settings, atoms))
