'''Python Fireworks API tests adapted from FireWorks'''
from fireworks import Firework, FWorker, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket
from fw_tutorials.firetask.addition_task import AdditionTask
from pathlib import Path
from sys import exit

from hilde.settings import Settings

def print_message(message):
    print(message)

# set up the LaunchPad and reset it
print("*******Testing Local FireWorks Installation*******")
print("Getting LaunchPad from file ~/.fireworks/my_launchpad.yaml.")

launchpad = LaunchPad.from_file(str(Path.home()) + "/.fireworks/my_launchpad.yaml")
print("*******Resetting LaunchPad*******")
launchpad.reset('', require_password=False)

# create the Firework consisting of a custom "Addition" task
firework = Firework(AdditionTask(), spec={"input_array": [1, 2]}, name="AdditionTestLocal")

# store workflow and launch it locally
launchpad.add_wf(firework)
launch_rocket(launchpad, FWorker())

print("*******Unmodified FireWorks is working locally*******")

print("*******Testing the local version of HiLDe's FireWorks modifications*******")

from hilde.tasks.fireworks.general_py_task import generate_firework
from hilde.tasks.fireworks.fw_action_outs import fireworks_no_mods_gen_function
from hilde.fireworks_api_adapter.launchpad import LaunchPadHilde

launchpad = LaunchPadHilde.from_file(str(Path.home()) + "/.fireworks/my_launchpad.yaml")

# define four individual FireWorks used in the Workflow
fw = generate_firework(
    print_message,
    fireworks_no_mods_gen_function,
    func_kwargs={},
    fw_settings={"fw_name": "testing", "spec":{}},
    args=["\nTesting HiLDe FireWorks API\n"],
    inputs=[],
)

# store workflow and launch it locally
launchpad.add_wf(fw)
launch_rocket(launchpad, FWorker())
print("*******The local version of HiLDe's FireWorks modifications work*******")

print("*******Testing the remote connection with FireWorks*******")

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
settings = Settings()
Si.set_calculator(setup_aims(settings=settings))

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
    fireworks_no_mods_gen_function,
    func_kwargs={"workdir": wd},
    atoms = Si,
    calc = Si.calc,
    fw_settings=fw_settings,
)
fw2 = generate_firework(
    print_message,
    fireworks_no_mods_gen_function,
    func_kwargs={},
    fw_settings={"fw_name": "testing", "spec":{}},
    args=["\n Connection successful \n"],
    inputs=[],
)
workflow = Workflow([fw, fw2], {fw: [fw2]})
launchpad.add_wf(workflow)
fworker = FWorker()
controlpath = str(Path.home() / ".ssh/sockets" / (settings.fireworks.remote_user + "@$cobra.mpcdf.mpg.de-22"))
print(f"Using the ssh ControlPath: {controlpath}")
rapidfire(launchpad, fworker, strm_lvl="INFO", remote_host=['cobra01i.mpcdf.mpg.de'], reserve=True, gss_auth=True, controlpath=controlpath)

print("\n\n*******All Tests Successful*******\n\n")
