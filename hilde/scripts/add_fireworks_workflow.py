from ase.io import read
from argparse import ArgumentParser

from fireworks import Workflow

from hilde import DEFAULT_CONFIG_FILE
from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.settings import Settings
from hilde.workflows.workflow_generator import generate_workflow

def main():
    parser = ArgumentParser(description="create a configuration file and workdir")
    parser.add_argument("workflow", type=str, help="Workflow description file")
    args = parser.parse_args()

    workflow = Settings(settings_file=args.workflow)
    atoms = read(workflow.geometry.file)
    steps = []

    for step_file in workflow.workflow.step_files:
        steps.append(Settings(settings_file=step_file))

    launchpad = LaunchPadHilde.from_file("/home/purcell/.fireworks/my_launchpad.yaml")
    launchpad.add_wf(generate_workflow(steps, atoms))
