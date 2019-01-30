"""Tests the ability of FireWorks to connect to a remote host"""
# Import modules from ase and fireworks
from pathlib import Path

from ase.build import bulk

from fireworks import FWorker, Workflow

from hilde.fireworks.combined_launcher import rapidfire
from hilde.fireworks.launchpad import LaunchPadHilde
from hilde.fireworks.workflow_generator import generate_firework
from hilde.settings import Settings
from hilde.templates.aims import setup_aims
from hilde.tasks.calculate import calculate as hilde_calc
from hilde.tasks.fireworks.fw_out.general import fireworks_no_mods_gen_function
from hilde.tasks.fireworks.general_py_task import TaskSpec


def print_message(message):
    """Prints a message"""
    print(message)


print("*******Testing the remote connection with FireWorks*******")

# Intialize Structures and database
Si = bulk("Si", "diamond")
settings = Settings()
Si.set_calculator(setup_aims(settings=settings))

# Port changes are for my setup
launchpad = LaunchPadHilde.auto_load()
fw_settings = {
    "fw_name": "test_remote",
    "spec": {
        # Submission script changes are controled by the _queueadapter dictionary
        "_queueadapter": {
            # Keys are the same that you define in "my_qadapter.yaml"
            "walltime": "00:01:00",
            "nodes": 1,
            "queue": "express",
            "ntasks_per_node": 32,
        }
    },
}

fw = generate_firework(
    atoms=Si,
    calc=Si.calc,
    fw_settings=fw_settings,
    func=hilde_calc,
    func_fw_out=fireworks_no_mods_gen_function,
    func_kwargs={},
)
task_spec = TaskSpec(
    print_message,
    fireworks_no_mods_gen_function,
    None,
    {},
    {},
    args=["\n Connection successful \n"],
)
fw2 = generate_firework(task_spec, fw_settings={"fw_name": "testing", "spec": {}})
workflow = Workflow([fw, fw2], {fw: [fw2]})
launchpad.add_wf(workflow)
fworker = FWorker()
rapidfire(launchpad, fworker, strm_lvl="INFO", reserve=True, gss_auth=True)

print("\n\n*******All Tests Successful*******\n\n")
