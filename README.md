**Installation**

`python setup.py install --prefix /your/preferred/folder`

**Basic Setup**

`cp hilde.cfg.template hilde.cfg`

and edit according to system.

**New Features**
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

**Setup of Docker Images to Speed Up Testing**
```
mkdir -p docker
cp Dockerfile requirements.txt docker
cd docker
docker login registry.gitlab.com
docker build -t registry.gitlab.com/flokno/hilde .
docker push registry.gitlab.com/flokno/hilde
```
**Setup of FireWorks on Computational Resources**
* FireWorks on the clusters
  * Download/clone from https://github.com/materialsproject/fireworks.git and move to that directory
  * Modify fw\_tuotirals/woker/my\_fworker.yaml and fw\_tuotirals/woker/my\_launchpad.yaml and
    copy them to $HOME/.fireworks
    * If your connection to the MongoDB is not secure (limited IP/port Access) set up a user/password for it and place those credentials in my_launchpad.yaml
    * Do not use your cluster/computer credentials for this
  * Modify the correctfw\_tuotirals/queue\_???.yaml file for your submission system
    and copy it to $HOME/.fireworks/my\_qadapter.yaml
  * Find the FireWorks install directory with lpad version and modify
    $FW_INSTALL_DIR/fireworks/fw_config.py:
    * LAUNCHPAD_LOC: $HOME/.fireworks/my_launchpad.yaml
    * FWORKER_LOC: $HOME/.fireworks/my_fworker.yaml
    * QUEUEADAPTER_LOC:$HOME/.fireworks/my_qadapter.yaml
  * Keep all PyTask functions up to date on all machines, there are no consistency checks
* Setup a MongoDB database for fireworks
  * Best to have it always accessible by all machines that need it
  * Check with the cluster management on what solution they'd prefer
  * Modify all my\_launchpad.yaml files such that they point to the database
* Connections between computers
  * Passwordless connections are preferred
  * If this is not possible you can pass the password as a command line argument, (delete history afterwards)
* FireWorks Etiquette
  * Name all Fireworks/WorkFlows
  * If you are using a shared launchpad only use lpad reset if everyone is okay with that
  * Keep the launch pad clean, once Workflows and FireWorks are complete/you don't need them,
  delete them
