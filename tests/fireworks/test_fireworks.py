""" An example of how to use FireWorks in conjunction with HilDe"""
import shutil
from pathlib import Path

from ase.build import bulk
from vibes.context import TaskContext
from vibes.fireworks.launchpad import LaunchPad
from vibes.fireworks.rocket_launcher import rapidfire
from vibes.fireworks.workflows.workflow_generator import generate_workflow
from vibes.helpers.paths import cwd
from vibes.settings import Settings


parent = Path(__file__).parent


def test_fireworks():
    settings = Settings(settings_file=parent / "workflow.in")
    settings.calculator["mkdir"] = False
    wflow = TaskContext(name=None, settings=settings)
    atoms = bulk("Ni", "fcc", a=3.5)
    wflow.atoms = atoms

    atoms.set_calculator(wflow.calculator)

    lp = LaunchPad(strm_lvl="INFO")
    lp.reset("", require_password=False)

    workflow = generate_workflow(wflow, atoms, None, True)
    lp.add_wf(workflow)

    with cwd("fireworks_launchers", mkdir=True):
        rapidfire(
            lp,
            wflow_id=workflow.root_fw_ids,
            strm_lvl="INFO",
            sleep_time=2,
            max_loops=2,
        )
    shutil.rmtree("test_run/")
    shutil.rmtree("fireworks_launchers/")

    assert lp.get_wf_by_fw_id(1).state == "COMPLETED"


if __name__ == "__main__":
    test_fireworks()
