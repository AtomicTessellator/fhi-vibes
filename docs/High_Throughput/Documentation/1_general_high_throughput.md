```
[fireworks]
name = example_run
config_dir: "~/.fireworks"
tasks2queue = ["vibes.relaxation.bfgs.relax", "vibes.fireworks.tasks.calculate_wrapper.wrap_calc_socket", "vibes.k_grid.converge_kgrid.converge_kgrid", "vibes.fireworks.tasks.calculate_wrapper.wrap_calculate", "vibes.fireworks.tasks.md.run" ]

[fireworks.workdir]
remote = "test_run/"
local   = "test_run/"

[fireworks.remote]
host = ["remote.host.path"]
config_dir = ["/path/to/remote/home/.fireworks/"]
launch_dir = "."

[fireworks.remote.authorization]
user = remote_username
password = null

[fireworks.remote.launch]
njobs_queue = 0
njobs_block = 500
reserve = True
nlaunches = 0
sleep_time = 60
```

## Sections

### `[fireworks]`
General parameters for the FireWorks workflows

#### `name`
`str`: The name that will be perpended the workflow to better organize the LaunchPad

#### `config_dir`
`str`: Directory where FireWorks configuration file are located (Default set in `.fireworksrc` file)

#### `tasks2queue`
`list(str)`: List of functions to send to the queue (Default set in `.fireworksrc` file)

### `[fireworks.workdir]`

These are used to define the base working directory on remote and local machines

#### `local`

`str`: Base working directory on your local machine

#### `remote`

`str`: Base working directory on a remote directory (Default is `fireworks.workdir.local`)

### `[fireworks.remote]`

Parameters for remote FireWorks workers (Default defined in `.fireworksrc`)

#### `host`

`list(str)`: List of remote hosts to send jobs to (Default defined in `.fireworksrc`)

#### `config_dir`

`list(str)`: List of remote FireWorks configuration directories (Default defined in `.fireworksrc`)

#### `launch_dir`

`str`: Default launch directory on the remote host (Default defined in `.fireworksrc`)

### `[fireworks.remote.authorization]`

Parameters for authentication remote FireWorks workers (Default defined in `.fireworksrc`)

#### `user`

`str`: remote host username (Default defined in `.fireworksrc`)

#### `password`

`str`: remote host password (not recommended) (Default defined in `.fireworksrc`)

### `[fireworks.remote.launch]`

Parameters for launching jobs on FireWorks workers (Default defined in `.fireworksrc`)

#### `njobs_queue`

`int`: number of jobs to have on the queue at any given time (0 no limit) (Default defined in `.fireworksrc`)

#### `njobs_block`

`int`: number of launches to have in a single FireWorks block directory (Default defined in `.fireworksrc`)

#### `reserve`

`bool`:  If True run FireWorks in reservation mode (Default defined in `.fireworksrc`)

#### `nlaunches`

`int`: Maximum number of jobs to launch at any given (0 no limit) (Default defined in `.fireworksrc`)

#### `sleep_time`

`float`: Time to sleep in seconds between checking for jobs to run (Default defined in `.fireworksrc`)
