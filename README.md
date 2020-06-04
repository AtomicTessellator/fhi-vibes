FHI-vibes
===

Welcome to `FHI-vibes`, a `python` package for _ab initio_ modeling of vibrational properties in anharmonic solids.

## Overview

- [Documentation](Documentation/0_intro/)
- [Tutorial](Tutorial/0_intro/)
- If you are interested in scientific work that was performed using `FHI-vibes`, please have a look at [References](References.md)

`FHI-vibes` is preparing a submission to [JOSS](https://joss.theoj.org/).

## Installation

### Prerequisites

- A working `python3.7+` or `python3.6` (see remarks below) environment, e.g., provided by [anaconda](https://docs.conda.io/en/latest/miniconda.html)

- A working `fortran` compiler, e.g., obtained by
  
  - `apt-get install gfortran` in Debian-derived systems, or
  - `conda install -c conda-forge fortran-compiler` when `conda` is used.


### Install `vibes`

`FHI-vibes` can be installed simply via pip:

```bash
pip install fhi-vibes
```

**(Important: If you run in to version conflicts that you cannot solve, use a virtual environment created with `python -m venv` or `conda create`.)**

### Configuration

Configure `vibes` by creating a `~/.vibesrc` configuration file in the home directory. To this end, first run

```
vibes template configuration > ~/.vibesrc
```

and edit according to system. The `aims_command` is a command or script that takes care of running aims. This can be either just `mpirun aims.x`, or a script loading necessary modules etc. and finally calling `srun aims.x` on a cluster.

**You're now good to go!** Just make sure your vibes virtual environment is activated.

### Remarks for `python3.6`

On `python3.6`, please install `importlib_resources` and `dataclasses` via 

```bash
pip install importlib_resources dataclasses
```

### Autocompletion

To activate autocompletion of `vibes` subcommands, add this to your `.bashrc`:

```bash
eval "$(_VIBES_COMPLETE=source vibes)"
```

and source it.

If you use the `fishshell`, add a file `~/.config/fish/completions/vibes.fish` containing

```bash
eval (env _VIBES_COMPLETE=source-fish vibes)
```
