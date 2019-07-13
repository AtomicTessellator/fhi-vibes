hilde
===
## Installation

External dependencies:

```bash
apt-get install gfortran liblapack-dev liblapacke-dev mongodb
```

Make sure `poetry` is installed:
```
pip install poetry
```

Create a virtual environment
```bash
python3 -m venv hilde_venv
source hilde_venv/bin/activate
```

Install `hilde` with `poetry`
```bash
poetry install
```

Configure Hilde by creating a `~/.hilderc` configuration file in the home directory. To this end, first

```
hilde template configuration
```

and edit according to system. The `aims_command` is a command or script that takes care
of running aims. This can be either just `mpirun aims.x`, or a script loading necessary
modules etc. and finally calling `srun aims.x` on a cluster.

Then copy to your home folder with 
```
cp hilderc ~/.hilderc
```

**You're now good to go!** Just make sure your hilde virtual environment is activated.

### Autocompletion
To activate autocompletion of `hilde` subcommands, add this to your `.bashrc`:
```bash
eval "$(_HILDE_COMPLETE=source hilde)"
```
and source it.

If you use the `fishshell`, add a file `~/.config/fish/completions/hilde.fish` containing
```bash
eval (env _HILDE_COMPLETE=source-fish hilde)
```

### Remarks for Cluster
On clusters with `conda` environment, it is typically best to have a minimal base
environment that holds nothing but `python`, `numpy`, and `scipy` to benefit from `mkl`
speedup:
```bash
conda create -n py37 python=3.7 numpy scipy mkl
```
From within the conda environment (`conda activate py37`), a `venv` can be created that
holds the other dependencies. To benefit from `numpy` and `mkl`, the `venv` can be
created with
```bash
python3 -m venv hilde_venv --system-site-packages
source hilde_venv/bin/activate
```
One should verify that the `numpy` and `scipy` packages are indeed the ones from the
conda environment by inspecting 
```bash
conda list | grep numpy
```
 within conda and
```bash
pip list | grep numpy
```
with `hilde_venv`. To enforce this, the `pyproject.toml` file
can be adjusted.

Don't forget to activate `hilde_venv` before starting to work!

## Settings Files

`hilde` uses the Python `configparser` module for parsing settings files named
`settings.in` and the configuration file `.hilderc`. The
parser is augmented by `JSON` so it understands any input on the right hand side that is
valid `JSON`. The inputs get converted to Python objects according to [this conversion
table](https://realpython.com/python-json/#serializing-json).

**New Features**
* Simplified Settings Files:
  * Settings files named `settings.in` are automatically parsed when calling
    `Settings()` within Hilde.
  * The configuration file `hilde.cfg` gets installed to system.
* Molecular dynamics workflow with input and output files
  * see hilde/examples/md
* Phonopy workflow with input and output files
  * see hilde/examples/phonopy
* Relaxation workflow with config file and output files
  * see hilde/examples/relaxation
* YAML Trajectories:
  * save MD trajectories as YAML with tools in `hilde.trajectories`
  * example in `hilde/examples/trajectory/trajectory.son.ipynb`
* Emails:
  * send notifications via email with `hilde.helpers.notifications.send_simple_mail`
* Watchdogs:
  * supervise e.g. an MD to estimate when the walltime will be reached.
    Example in `examples/md/md_with_watchdog.ipynb`
* Wrapper for `phono3py`
  * Preprocess and re-creation of Phono3py objects from precomputed force
  constants, see examples
* Wrapper for `phonopy`
  * Preprocess and (some) postprocess, see examples
* Templates
  * `from hilde.templates.lammps import setup_lammps_si` to provide lammps calculator
* Brillouin zone helpers
  * `hilde.helpers.brillouinzone` features `get_paths`, `get_bands`, and
  `get_labels` to provide paths in the BZ that can be fed to `phonopy` via
  `phonon.set_bandstructure(bands)`, and
  `phonon.plot_band_structure(labels=labels)`.
  * These functions are used by `hilde.phonopy.plot_dos_and_bandstructure` to
  plot DOS and bandstructure in the working directory.
* Scripts:
  * `make_supercell`: create supercell from supercell matrix or
  target target
  * `geometry_info`: print geometry information for given input
  structure
* Symmetry Block Generation Functions
  * `AtomsInput`: A storage class that stores relevant information about a structure
  * `write_sym_constraints_geo`: Read any geometry.in file and use the list of `AtomInputs`
  to create a new supercell with a user defined symmetry block added to it
* FireWorks integration
  * Functions that can be used with PyTask to use FireWorks as a job manager
  * Jobs can now be submitted to the queue from a local machine and have the results processed locally


**Setup of FireWorks on Computational Resources**

See also: `doc/README_FHI_FireWorksConnections.md`
* Overview of Managing FireWorks Remotely
  * FireWorks does not copy functions but only finds them in the PYTHONPATH
  * To pass it functions give it the function_module.function_name as a str
  * Functions that are run on the local machine
    * All functions/files that set up FireWorks
      * All scripts that initially call hilde.fireworks.tasks.generate_firework
      * .cfg Files that define the steps (if used)
      * All functions used by a Fireworks without a task that calls a function in task2queue list
    * claunch_hilde and associated functions
  * Function that are run on the remote machine
    * All functions used by a Firework with a task that calls a function in task2queue
    * qluanch_hilde and associated functions
  * Functions that can run on both machines
    * All FireWorks API functions
    * All database accessors functions
    * Spec modifying functions (hilde.fireworks.tasks.fw_action_outs)
    * hilde.fireworks.tasks.generate_firework
  * Machine specific settings such as the aims_command is handled dynamically
    * It automatically changes when called on a machine
    * Can always use local settings without an issue
* Prerequisites for using FireWorks
  * Fabric 2 (for remote connections)
  * paramiko (used by Fabric 2)
  * python-gssapi (for gss authorization)
  * pymongo
* Using FireWorks on the clusters
  * Download/clone from https://github.com/materialsproject/fireworks.git and move the FireWorks directory
  * Modify fw\_tutorials/worker/my\_fworker.yaml and copy it to $HOME/.fireworks
    * Probably do not need to do any modifications if running on similar environments
    * Useful if you want to run specific jobs on specific machines without specified reservations
  * Modify fw\_tutorials/worker/my\_launchpad.yaml and copy it to $HOME/.fireworks
    * host: Host to the DB server
      * If connected through an ssh tunnel use localhost
    * port: Port the DB server is listening on
      * If connected through an ssh tunnel use the port connected the DB server via the tunnel
    * username: username used to access the database
    * password: password used to access the database
    * logdir: default directory to store logs
    * strm_lvl: How much information the launchpad prints by default
  * Modify the correct fw\_tutorials/queue\_???.yaml file for your submission system and copy it to $HOME/.fireworks/my\_qadapter.yaml
    * Only used on clusters
    * Set to minimal queue defaults
      * nodes = 1
      * ntasks_per_node = 32
      * walltime = "00:30:00"
      * queue = "express"
      * logdir = /some/path/that/must/exist (make sure this exists)
  * Find the FireWorks install directory with lpad version and modify
    $FW_INSTALL_DIR/fireworks/fw_config.py:
    * LAUNCHPAD_LOC: $HOME/.fireworks/my_launchpad.yaml
    * FWORKER_LOC: $HOME/.fireworks/my_fworker.yaml
    * QUEUEADAPTER_LOC: $HOME/.fireworks/my_qadapter.yaml
* Setup a MongoDB database for fireworks
  * Best to have it always accessible by all machines that need it
  * Check with the cluster management on what solution they'd prefer
* Connections between computers
  * Passwordless connections are preferred
  * If this is not possible you can pass the password as a command line argument, (delete
    bash history afterwards)
* FireWorks Etiquette
  * Name all Fireworks/WorkFlows
  * If you are using a shared launchpad only use lpad reset if everyone is okay with that
