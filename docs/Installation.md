# Installation

## Prerequisites

- A working`python3.6+` environment, e.g., provided by [anaconda](https://docs.conda.io/en/latest/miniconda.html).

- A working `fortran` compiler, e.g., obtained by:
    - `apt-get install gfortran` in Debian-derived systems, or
    - `conda install -c conda-forge fortran-compiler` when `conda` is used.

- If you want to use `FHI-aims` for running _ab initio_ calculations, make sure you have a recent version that supports the iPi socket communication (this is the default for any version newer than the `200112_2` release when using the [CMake build system](https://aims-git.rz-berlin.mpg.de/aims/FHIaims/-/wikis/CMake-Tutorial)).


## Install `vibes`

`FHI-vibes` can be installed simply via pip:

```bash
pip install --user fhi-vibes
```

The `--user` option makes sure that the installation occurs in your homefolder under `~/.local/bin`, as typically
necessary on computing clusters. Please make sure that the `~/.local/bin` folder is listed in your `PATH`.

**If you run into problems, please have a look at our [troubleshooting section.](#Troubleshooting)**

## Configuration

Configure `vibes` by creating a `~/.vibesrc` configuration file template in the home directory. To this end, first run

```
vibes template configuration vibes > ~/.vibesrc
```

and edit the configuration file as described below:

### `basissetloc`

The `basissetloc` should point to the folder containing FHI-aims' species defaults, e.g., the `/path/to/FHIaims/species_defaults` folder.

### `aims_command`

The `aims_command` should be an executable script that takes care of setting up the environment and then running FHI-aims, for example a file called `/path/to/FHIaims/run_aims.sh`  that looks roughly like this (depends on you system!):

```
#!/bin/bash -l

ulimit -s unlimited
export OMP_NUM_THREADS=1

module purge
module load intel impi mkl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MKLROOT}/lib/intel64/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${INTEL_HOME}/lib/intel64/

srun /path/to/FHIaims/build/aims.x
```

The script `run_aims.sh` has to be marked as executable, e.g., by running `chmod +x /path/to/FHIaims/run_aims.sh`.

**You're now good to go!**

## Autocompletion

To activate autocompletion of `vibes` subcommands, add this to your `.bashrc`:

```bash
eval "$(_VIBES_COMPLETE=source vibes)"
```

and source it.

If you use the `fishshell`, add a file `~/.config/fish/completions/vibes.fish` containing

```bash
eval (env _VIBES_COMPLETE=source-fish vibes)
```



## Troubleshooting

- `ModuleNotFoundError: No module named 'importlib_resources'`
    - Solution: `pip install importlib_resources dataclasses`
- `RuntimeError: Click will abort further execution because Python 3 was configured to use ASCII as encoding for the environment. Consult https://click.palletsprojects.com/python3/ for mitigation steps`
    - Solution:  `export LC_ALL=C.UTF-8 ; export LANG=C.UTF-8`
- `-bash: vibes: command not found`
    - Solution: `export PATH=$PATH:~/.local/bin`
- `ImportError: numpy.core.multiarray failed to import`
    - Solution: `pip install numpy -U` (or `conda update numpy` if you use conda)
- Various version conflicts
    - Consider using a [virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), e.g., via `conda create -n py38 -c anaconda python=3.8 numpy scipy`

- `/tmp/pip-build-env-xak_2vfy/overlay/lib/python3.7/site-packages/numpy/core/_multiarray_umath.cpython-37m-x86_64-linux-gnu.so: failed to map segment from shared object: Operation not permitted`
    - This might happen on HPC systems with limited access rights. The solution is to provide a writable `tmp` folder, e.g. via `mkdir ~/condatmp && export TMPDIR=~/condatmp/`

If your problem is not listed here, please [file an issue in our issue tracker](https://gitlab.com/vibes-developers/vibes/-/issues).
