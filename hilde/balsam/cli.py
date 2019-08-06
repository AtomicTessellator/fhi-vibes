"""click interface for the balsam"""
# Importing to ensure balsam is installed
import balsam.launcher.dag as dag

import click

from ase.io.aims import read_aims

from pathlib import Path

from hilde.balsam.apps import register
from hilde.balsam.workflow.workflow_generator import generate_workflow

from hilde.aims.context import AimsContext
from hilde.aims.setup import setup_aims
from hilde.cli.misc import AliasedGroup
from hilde.settings import TaskSettings, AttributeDict, Settings


@click.command(cls=AliasedGroup)
def balsam():
    """Access HiLDe's balsam wrappers"""
    pass


@balsam.command("add_wf")
@click.option("-w", "--workflow", default="workflow.in")
def add_wf(workflow):
    """Adds a workflow to the launchpad"""
    register()
    settings = Settings(settings_file=workflow, read_config=False)
    wflow = TaskSettings(name=None, settings=settings, read_config=False)

    structures = []

    if "geometry" in wflow:
        if "files" in wflow.geometry:
            for file in Path.cwd().glob(wflow.geometry.files):
                structures.append(read_aims(file))
        if "file" in wflow.geometry:
            structures.append(wflow.get_atoms())
    else:
        raise IOError("No geometry file was specified")

    settings_full = Settings(settings_file=workflow, read_config=True)
    if "basisset" not in settings_full:
        settings_full["basisset"] = AttributeDict({"type": "light"})
    calc = setup_aims(AimsContext(settings=settings_full))

    for atoms in structures:
        atoms.set_calculator(calc)
        settings = Settings(settings_file=workflow, read_config=False)
        wflow = TaskSettings(name=None, settings=settings, read_config=False)
        wflow.atoms = atoms

        if "basisset" not in wflow and "basisset" in wflow.general:
            wflow["basisset"] = AttributeDict({"type": wflow.general.basisset})
        elif "basisset" not in wflow:
            wflow["basisset"] = AttributeDict({"type": "light"})

        generate_workflow(wflow)
