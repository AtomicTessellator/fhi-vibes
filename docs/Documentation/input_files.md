# Input File Format

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

#### Example

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
    "calculator": {"name": "lj", "parameters": {"sigma": 3.4}},
    "files": {"geometry": "geometry.in"},
    "relaxation": {
        "driver": "BFGS",
        "fmax": 0.001,
        "kwargs": {"maxstep": 0.2},
        "workdir": "lj.relaxation",
    },
}
```

### `[files]` Section

This section contains filenames. 

#### `geometry`

`geometry` gives the name of the geometry input file to be used for a calculation:

```python
file = settings.files.get("geometry")

atoms = ase.io.read(file)
```

If there is just one geometry necessary for the given task and it is stored in `geometry.in`, this section can be omitted altogether.

#### `geometries`

Via `geometries`, a wildcard expression for finding geometry files for computation can be given, e.g. `geometries: samples/geometry.in.*` would specifiy to run a calculation for all geometry input files found in the folder `samples`.

```python
files = sorted(glob(settings.files.get("geometries")))

atoms_to_compute = [ase.io.read(file) for file in files]
```


#### `primitive`
Give a reference primitive cell in a file, e.g., `primitive: geometry.in.primitive`

#### `supercell`
Give a reference supercell in a file, e.g., `supercell: geometry.in.supercell`

#### Example
Example for specifying to run a job for the structure in `geometry.in`, while attaching a reference primitive and supercell to the output trajectory:

```
[files]
geometry:                      geometry.in
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell

...
```