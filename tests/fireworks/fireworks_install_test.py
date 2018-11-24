'''Python Fireworks API tests adapted from FireWorks'''
from fireworks import Firework, FWorker, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket
from fw_tutorials.firetask.addition_task import AdditionTask
from pathlib import Path
from sys import exit

def print_message(message):
    print(message)

# set up the LaunchPad and reset it
print("Getting LaunchPad from file ~/.fireworks/my_launchpad.yaml.")

launchpad = LaunchPad.from_file(str(Path.home()) + "/.fireworks/my_launchpad.yaml")
launchpad.reset('', require_password=False)

# create the Firework consisting of a custom "Addition" task
firework = Firework(AdditionTask(), spec={"input_array": [1, 2]})

# store workflow and launch it locally
launchpad.add_wf(firework)
launch_rocket(launchpad, FWorker())

print("The native FireWorks installation is working, retrying with modified versions")

print("Testing the local version of HiLDe's FireWorks modification")

from hilde.tasks.fireworks.general_py_task import generate_firework
from hilde.tasks.fireworks.fw_action_outs import return_null_general
from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde

launchpad = LaunchPadHilde.from_file(str(Path.home()) + "/.fireworks/my_launchpad.yaml")

# define four individual FireWorks used in the Workflow
fw = generate_firework(
    print_message,
    return_null_general,
    func_kwargs={},
    fw_settings={"fw_name": "testing", "spec":{}},
    args=["\nTesting HiLDe FireWorks API\n"],
    inputs=[],
)

# store workflow and launch it locally
launchpad.add_wf(fw)
launch_rocket(launchpad, FWorker())

print("Testing the remote connection with FireWorks")

# Import modules from ase and fireworks
import numpy as np

from ase.build import bulk
from ase.calculators.emt import EMT
from ase.db.core import connect

from fireworks import Firework, FWorker, PyTask, FWAction, Workflow

# Combined local/remote queue launching
from hilde.fireworks_api_adapter.combined_launcher import rapidfire
from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde
from hilde.templates.aims import setup_aims
# Import the hilde calculate function so both the local and remote machines have the same function in their path
from hilde.tasks.calculate import calculate as hilde_calc

# Intialize Structures and database
Si = bulk("Si", 'diamond')
Si.set_calculator(setup_aims())

# Port changes are for my setup
launchpad = LaunchPadHilde.from_file(str(Path.home()) + "/.fireworks/my_launchpad.yaml")
fw_settings = {
    "fw_name": "test_remote",
    "spec":
    {
        # Submission script changes are controled by the _queueadapter dictionary
        "_queueadapter":
        {
            # Keys are the same that you define in "my_qadapter.yaml"
            "walltime": "00:01:00",
            "nodes": 1,
        }
    },
}
wd = "/u/tpurcell/.fireworks/Si/"
fw = generate_firework(
    hilde_calc,
    return_null_general,
    func_kwargs={"workdir": wd},
    atoms = Si,
    calc = Si.calc,
    fw_settings=fw_settings,
)
fw2 = generate_firework(
    print_message,
    return_null_general,
    func_kwargs={},
    fw_settings={"fw_name": "testing", "spec":{}},
    args=["\n Connection successful \n"],
    inputs=[],
)
workflow = Workflow([fw, fw2], {fw: [fw2]})
launchpad.add_wf(workflow)
fworker = FWorker()

rapidfire(launchpad, fworker, strm_lvl="INFO", reserve=True, gss_auth=True)
