FHI-vibes
===

Welcome to `FHI-vibes`, a `python` package for _ab initio_ modeling of vibrational properties in anharmonic solids. `FHI-vibes` is intended to bridge between different methodologies, so to allow for a seamless assessment of vibrational properties with different approaches, ranging from the harmonic approximation to anharmonic MD.

## Overview

- [Tutorial](https://vibes-developers.gitlab.io/vibes/Tutorial/0_intro/)
- [Documentation](https://vibes-developers.gitlab.io/vibes/Documentation/0_intro/)
- If you are interested in scientific work that was performed using `FHI-vibes`, please have a look at [References](https://vibes-developers.gitlab.io/vibes/References/)

## Credits

`FHI-vibes` would not be possible without the following packages:

- The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
- [FHI-aims: FHI _ab initio_ molecular simulations](https://aimsclub.fhi-berlin.mpg.de/)
- [Phonopy](https://atztogo.github.io/phonopy/) and [Phono3py](https://atztogo.github.io/phono3py/)
- [hiPhive — High-order force constants for the masses](https://hiphive.materialsmodeling.org/index.html)

### How to cite these packages:

Please make sure to give credit to the right people when using `FHI-vibes`:

- [How to cite ASE](https://wiki.fysik.dtu.dk/ase/faq.html#how-should-i-cite-ase)
- [How to cite FHI-aims](https://aimsclub.fhi-berlin.mpg.de/aims_publications.php)
- [How to cite phonopy](https://phonopy.github.io/phonopy/citation.html)
- [How to cite phono3py](https://phonopy.github.io/phono3py/citation.html)
- [How to cite hiphive](https://hiphive.materialsmodeling.org/credits.html)
- [How to cite FHI-vibes: coming soon]()

## Installation

### Prerequisites

- A working `python3.7+` or `python3.6` (see remarks below) environment, e.g., provided by [anaconda](https://docs.conda.io/en/latest/miniconda.html)

- A working `fortran` compiler, e.g., obtained by
  
  - `apt-get install gfortran` in Debian-derived systems, or
  - `conda install -c conda-forge fortran-compiler` when `conda` is used.
- If you want to use `FHI-aims` for running _ab initio_ calculations, make sure you have a recent version that supports the iPi socket communication.


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
