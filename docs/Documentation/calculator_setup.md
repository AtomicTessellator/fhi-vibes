Calculator Setup
===

## Example

```
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4
```

## Sections

### `[files]`

This section contains filenames. 

#### `geometry`

`geometry` gives the name of the geometry input file to be used for a calculation:

```python
file = settings.files.get("geometry")

atoms = ase.io.read(file)
```



### `[calculator]`

This section specifies which `ase.Calculator` should be set up and how.

#### `name`

The name of the [ASE calculator class name](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#supported-calculators).

### `[calculator.parameters]`

These keywords are used 1:1 to set up the ASE calculator:

```python
cls = get_calculator_class(settings.calculator.get("name"))

calculator = cls(**settings.calculator.get("parameters"))
```

## More Options for `FHI-aims`

```
[calculator]
name:                          aims

[calculator.parameters]
xc:                            pw-lda

[calculator.kpoints]
density:                       3.5

[calculator.basissets]
default:                       intermediate
fallback:                      light
# system specific
# O:                           tight
# Ga:                          intermediate

[calculator.socketio]
port:                          12345
```

### `[calculator.parameters]`

These keywords correspond one-to-one to the FHI-aims keywords that  are written to `control.in`.

### `[calculator.kpoints]`

#### `density`

Compute `k_grid` such that the density of kpoints does not fall below this value in $\require{mediawiki-texvc} \AA^{-3}$ .

### `[calculator.basissets]`

Specify which basissets to use.

#### `default`

The default basis set to use, can be `light`, `intermediate`,`tight`, or `really_tight`.

#### `fallback`

The fallback option in case the specified basis set could not be found (`intermediate` basis sets are currently not compiled for each element)

#### Species dependent

The basis set can be given per chemical species by including the species and its desired basis set (uncomment, e.g., `O` in the example above.)

### `[calculator.socketio]`

Set up socket communication via [`SocketIOCalculator`](https://wiki.fysik.dtu.dk/ase/ase/calculators/socketio/socketio.html?highlight=socketio#ase.calculators.socketio.SocketIOCalculator)

#### `port`

The socket port to use.

- `null`: don't use the socket.
- `0`-`65535`: use this port.

