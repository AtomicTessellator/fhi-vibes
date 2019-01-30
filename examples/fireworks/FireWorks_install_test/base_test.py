"""Tests the base FireWorks Installation"""

from fireworks import Firework, FWorker, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket
from fw_tutorials.firetask.addition_task import AdditionTask


# set up the LaunchPad and reset it
print("*******Testing Local FireWorks Installation*******")
print("Getting LaunchPad from file ~/.fireworks/my_launchpad.yaml.")

launchpad = LaunchPad.auto_load()
print("*******Resetting LaunchPad*******")
# launchpad.reset('', require_password=False)

# create the Firework consisting of a custom "Addition" task
firework = Firework(
    AdditionTask(), spec={"input_array": [1, 2]}, name="AdditionTestLocal"
)

# store workflow and launch it locally
launchpad.add_wf(firework)
launch_rocket(launchpad, FWorker())
print("*******Unmodified FireWorks is working locally*******")
