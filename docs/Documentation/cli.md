# Command Line Interface (CLI)

FHI-vibes comes with a command line interface (CLI) for

- creating input files (`vibes  template`),
- informing about input and output files (`vibes info`),
- running calculations (`vibes run` and `vibes submit`),
- processing calculations (`vibes output`), and
- performing several other tasks like converting output files (`vibes utils`).


## `vibes template`

Create template files

```
$ vibes template --help

Usage: vibes template [OPTIONS] COMMAND [ARGS]...

  provide template input files for tasks and workflows

Options:
  -h, --help         Show this message and exit.

Commands:
  calculator     Calculator templates: aims, lj
  configuration  Configuration templates: .vibesrc, .fireworksrc
  md             provide template input for MD simulation (default: NVE)
  phonopy        provide template input for phonopy workflow.
  relaxation     provide template input for relaxation workflow.
  slurm          provide template slurm settings
```

Each of the the sub-commands has its own `--help` for additional information.
The templates are printed to screen and can be piped to a file with `| tee`, `>` or `>>`.

## `vibes info`
```
$ vibes info --help

Usage: vibes info [OPTIONS] COMMAND [ARGS]...

  inform about content of a file

Options:
  -h, --help  Show this message and exit.

Commands:
  csv             show contents of csv FILE
  geometry        inform about a structure in a geometry input file
  md              inform about content of a settings.in file
  netcdf          show contents of netCDF FILE
  phonopy         inform about a phonopy calculation
  relaxation      inform about geometry optimization
  settings        inform about content of a settings file
  trajectory      inform about content of trajectory file
```