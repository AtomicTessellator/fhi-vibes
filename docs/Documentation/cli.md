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
  --allow_overwrite  [default: False]
  -h, --help         Show this message and exit.

Commands:
  calculator     Calculator templates: aims, lj
  configuration  Configuration templates: .vibesrc, .fireworksrc
  md             provide template md.in for molecular dynamics workflow.
  phonopy        provide template phonopy.in for phonopy workflow.
  relaxation     provide template relaxation.in for relaxation workflow.
  slurm          provide template slurm settings
```

