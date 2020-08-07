# Configuring FHI-vibes for interaction with FireWorks
Once FireWorks is setup and running  a few additional steps are needed to configure FHI-vibes to work with FireWorks.
If you are having trouble installing FireWorks see [their documentation](https://materialsproject.github.io/fireworks/) or [our installation guide](../Installation/0_setup.md).

## Creating a `.fireworksrc` file
The first step in configuring FHI-vibes for use with FireWorks.
This file sets up default values for running FireWorks utilities that will be consistent throughout all calculation.
To make this file run:

`vibes template configuration fireworks > ~/.fireworksrc`

```
[fireworks]
config_dir: "~/.fireworks" # Directory containing the *yaml files for FireWorks
tasks2queue = ["vibes.relaxation.bfgs.relax", "vibes.fireworks.tasks.calculate_wrapper.wrap_calc_socket", "vibes.k_grid.converge_kgrid.converge_kgrid", "vibes.fireworks.tasks.calculate_wrapper.wrap_calculate", "vibes.fireworks.tasks.md.run" ] # DO NOT CHANGE

[fireworks.remote]
host = ["remote.host.path"] # List of remote host names
config_dir = ["/path/to/remote/home/.fireworks/"] # List of remote FireWorks configuration directories
launch_dir = "." # Default launch directory on the remote host

[fireworks.remote.authorization]
user = remote_username # remote host username
password = null # remote host password (not recommended try to use password-less login)

[fireworks.remote.launch]
njobs_queue = 0 # Number of jobs to have on the queue at any given time (0 no limit)
njobs_block = 500 # Number of launches to have in a single FireWorks block directory
reserve = True # If True run FireWorks in reservation mode
nlaunches = 0 # Maximum number of jobs to launch at any given (0 no limit)
sleep_time = 60 # Time to sleep in seconds between checking for jobs to run
```
For a complete description of each of these parameters see the [full documentation](../../Documentation/1_general_high_throughput).

## Testing if it works
Now that your FireWorks installation should be working properly test it with the vibes FireWorks test in `test/fireworks/test_fireworks.py`.
If the test runs successfully then when you run `lpad get_wflows` you should get the following output (created on should match today's date/time):
```
{
    "state": "COMPLETED",
    "name": "Ni_6d2a2be5a5c1c4549639b55c5403b438b3b0ccf7--1",
    "created_on": "2020-03-13T11:55:30.357000",
    "states_list": "C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C"
}

```
If you see that you have successfully set up the high-throughput portions of vibes. To use this on clusters you need to repeat the steps in Installing/Testing FireWorks on each machine you plan to use it on.
