Input Files
===

## The `jconfigparser` syntax

`FHI-vibes` uses [`jconfigparser`](https://pypi.org/project/jconfigparser/) to parse input files. `jconfigparser` is an extension to the `python` [standard library `configparser`](https://docs.python.org/3/library/configparser.html) with the following additions:

- Nested section names separated with `.` ,
- values are parsed by `json`,
- repeated keywords possible.

## Example

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