"""Tests the HiLDe modifications to FireWorks"""
from pathlib import Path

from fireworks import FWorker
from fireworks.core.rocket_launcher import launch_rocket

from hilde.fireworks.workflow_generator import generate_firework
from hilde.fireworks.launchpad import LaunchPadHilde
from hilde.tasks.fireworks.fw_out.general import fireworks_no_mods_gen_function
from hilde.tasks.fireworks.general_py_task import TaskSpec


def print_message(message):
    """Prints a message"""
    print(message)


print("*******Testing the local version of HiLDe's FireWorks modifications*******")
launchpad = LaunchPadHilde.auto_load()

# define four individual FireWorks used in the Workflow
task_spec = TaskSpec(
    print_message,
    fireworks_no_mods_gen_function,
    None,
    {},
    {},
    args=["\nTesting HiLDe FireWorks API\n"],
)
fw = generate_firework(task_spec, fw_settings={"fw_name": "testing", "spec": {}})

# store workflow and launch it locally
launchpad.add_wf(fw)
launch_rocket(launchpad, FWorker())
print("*******The local version of HiLDe's FireWorks modifications work**********")
