FHI-vibes
===

Welcome to `FHI-vibes`, a `python` package for _ab initio_ modeling of vibrational properties in anharmonic solids.

## Installation

### Prerequisites

- A working `python3.7+` or `python3.6` (see remarks below) environment, e.g., provided by [anaconda](https://docs.conda.io/en/latest/miniconda.html)

- A working `fortran` compiler, e.g., obtained by
  
  - `apt-get install gfortran` in Debian-derived systems, or
  - `conda install -c conda-forge fortran-compiler` when `conda` is used.

- `poetry` needs to be installed
  
  - `pip install poetry`
  - Make sure that `poetry` will no creat a virtual environment of its own by
    - `poetry config virtualenvs.create false`
    - or create a [virtual environment]([https://docs.python.org/3/library/venv.html](https://docs.python.org/3/library/venv.html)

### Install `vibes`

Go to the `vibes` folder and install `vibes` with `poetry`

```bash
cd /path/to/vibes
poetry install
```

**(Important: If you run in to version conflicts that you cannot solve, use a virtual environment created with `python -m venv` or `conda create`.)**

### Configuration

Configure Hilde by creating a `~/.vibesrc` configuration file in the home directory. To this end, first run

```
vibes template configuration
```

and edit according to system. The `aims_command` is a command or script that takes care of running aims. This can be either just `mpirun aims.x`, or a script loading necessary modules etc. and finally calling `srun aims.x` on a cluster.

Then copy this file to your home folder with 

```
cp vibesrc ~/.vibesrc
```

**You're now good to go!** Just make sure your vibes virtual environment is activated.

## Remarks for `python3.6`

On `python3.6`, please install `importlib_resources` via 

```bash
pip install importlib_resources
```

### Autocompletion

To activate autocompletion of `vibes` subcommands, add this to your `.bashrc`:

```bash
eval "$(_HILDE_COMPLETE=source vibes)"
```

and source it.

If you use the `fishshell`, add a file `~/.config/fish/completions/vibes.fish` containing

```bash
eval (env _HILDE_COMPLETE=source-fish vibes)
```

### Remarks for Cluster

First make sure that `poetry` doesn't create virtual environments on its own:

```bash
poetry config virtualenvs.create false 
```

On clusters with `conda` environment, it is typically best to have a minimal base
environment that holds nothing but `python`, `numpy`, and `scipy` to benefit from `mkl`
speedup:

```bash
conda create -n py37 python=3.7 numpy scipy mkl
```

From within the conda environment (`conda activate py37`), you can install just like above. To enforce using the `mkl` `numpy`, you can find out the version with

```bash
conda list | grep numpy
```

 or

```bash
pip list | grep numpy
```

and adjust the `pyproject.toml` file to match the exact version number.

Don't forget to activate `py37` with `conda activate py37` before starting to work!

## Settings Files

`vibes` uses the Python `configparser` module for parsing settings files named
`settings.in` and the configuration file `.vibesrc`. The
parser is augmented by `JSON` so it understands any input on the right hand side that is
valid `JSON`. The inputs get converted to Python objects according to [this conversion
table](https://realpython.com/python-json/#serializing-json).