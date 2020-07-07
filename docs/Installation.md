# Installation

## Prerequisites

- A working `python3.7+` or `python3.6` environment, e.g., provided by [anaconda](https://docs.conda.io/en/latest/miniconda.html). 
	- On `python3.6`, please `pip install importlib_resources dataclasses`

- A working `fortran` compiler, e.g., obtained by
  
  - `apt-get install gfortran` in Debian-derived systems, or
  - `conda install -c conda-forge fortran-compiler` when `conda` is used.
- If you want to use `FHI-aims` for running _ab initio_ calculations, make sure you have a recent version that supports the iPi socket communication (this is the default when using the [CMake build system](https://aims-git.rz-berlin.mpg.de/aims/FHIaims/-/wikis/CMake-Tutorial)).


## Install `vibes`

`FHI-vibes` can be installed simply via pip:

```bash
pip install --user fhi-vibes
```

**Important: If you run in to version conflicts that you cannot solve, use a virtual environment created with `python -m venv` or `conda create`.**

## Configuration

Configure `vibes` by creating a `~/.vibesrc` configuration file in the home directory. To this end, first run

```
vibes template configuration vibes > ~/.vibesrc
```

and edit according to system.

### `basissetloc`

The `basissetloc` points to your `/path/to/FHIaims/species_defaults` folder.

### `aims_command`

The `aims_command` should be an executable script that  takes care of running aims, for example a file called `run_aims`  that looks roughly like this (depends on you system!):

```
#!/bin/bash -l

ulimit -s unlimited
export OMP_NUM_THREADS=1

module purge
module load intel impi mkl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MKLROOT}/lib/intel64/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${INTEL_HOME}/lib/intel64/

srun aims.x
```

This file has to live on your `PATH` and has to be made executable, e.g., via `chmod +x`. 

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
