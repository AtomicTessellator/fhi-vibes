from ase.io import read
from argparse import ArgumentParser

from fireworks import Workflow

from hilde import DEFAULT_CONFIG_FILE
from hilde.fireworks.launchpad import LaunchPadHilde
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.settings import Settings
from hilde.workflows.workflow_generator import generate_workflow

def main():
    parser = ArgumentParser(description="create a configuration file and workdir")
    parser.add_argument("workflow", type=str, help="Workflow description file")
    parser.add_argument("--make_abs_path", action="store_true", help="Make all paths absolute")
    args = parser.parse_args()

    workflow = Settings(settings_file=args.workflow)
    atoms = read(workflow.geometry.file)
    steps = []

    for step_file in workflow.workflow.step_files:
        steps.append(Settings(settings_file=step_file))

    fw_settings = {"to_launchpad": True}
    if "name" not in fw_settings:
        fw_settings["name"] = ""
    fw_settings["name"] += atoms.symbols.get_chemical_formula() + "_" + hash_atoms_and_calc(atoms)[0]
    generate_workflow(steps, atoms=atoms, fw_settings=fw_settings, make_abs_path=args.make_abs_path)
