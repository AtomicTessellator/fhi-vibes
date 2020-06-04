Input Files
===

## Geometry input files
`FHI-vibes` uses the `FHI-aims` geometry description format `geometry.in`. A detailed documentation of the file format can be found [here](https://doi.org/10.1016/j.cpc.2009.06.022).

### Example
An example `geometry.in` for fcc-silicon reads
```
lattice_vector 0.00 2.72 2.72
lattice_vector 2.72 0.00 2.72
lattice_vector 2.72 2.72 0.00

atom_frac 0.00 0.00 0.00 Si
atom_frac 0.25 0.25 0.25 Si
```

## Task input files
For performing a specific task, say, a geometry optimization, `FHI-vibes` employs single input files describing the task, e.g., a `relaxation.in` file describing the optimization according to the [documentation](relaxation.md).

### The `jconfigparser` syntax
The input files are parsed using [`jconfigparser`](https://pypi.org/project/jconfigparser/) . `jconfigparser` is an extension to the `python` [standard library `configparser`](https://docs.python.org/3/library/configparser.html) with the following additions:

- Nested section names separated with `.` ,
- values are parsed by `json`,
- repeated keywords possible,
- [`configparser.ExtendedInterpolation`](https://docs.python.org/3/library/configparser.html#configparser.ExtendedInterpolation) is used per default.

### Example

An example for an input file for running a geometry optimization:

```
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4

[relaxation]
driver:                        BFGS
fmax:                          0.001
workdir:                       ${calculator.parameters:xc}.relaxation

[relaxation.kwargs]
maxstep:                       0.2
```

This file will be parsed to a nested dictionary:

```
settings = {
    "files": {"geometry": "geometry.in"},
    "calculator": {
        "basissets": {"default": "light"},
        "name": "aims",
        "parameters": {"k_grid": [4, 4, 4], "xc": "pw-lda"},
    },
    "relaxation": {
        "driver": "BFGS",
        "fmax": 0.001,
        "kwargs": {"logfile": "relaxation.log", "maxstep": 0.2},
        "workdir": "pw-lda.relaxation",
    },
}
```