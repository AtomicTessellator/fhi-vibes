# Calculator Setup

FHI-vibes can set up any [ASE calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#module-ase.calculators) for performing a calculation by providing the  calculator class `name` and the respective `parameters` in the input file.

## Example

```
...
[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4
...
```

This would set up a [Lennard Jones calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/others.html#lennard-jones) with a `sigma` value of 3.4 and default parameters otherwise.

## Sections

### `[calculator]`

This section specifies which `ase.Calculator` should be set up and how.

#### `name`

The name of the [ASE calculator class name](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#supported-calculators).

### `[calculator.parameters]`

These are the keywords used 1:1 to set up the ASE calculator:

```python
cls = get_calculator_class(settings.calculator.get("name"))

calculator = cls(**settings.calculator.get("parameters"))
```

## Options for `FHI-aims`

FHI-vibes is most tightly integrated with the FHI-aims calculator and provides some extra features for performing _ab initio_ calculations with FHI-aims. A minimal input section to set up an FHI-aims calculator looks like this:

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

These keywords correspond one-to-one to the FHI-aims keywords that  are written to `control.in`. Keyword-only arguments like `vdw_correction_hirshfeld` or `use_gpu` should be given with the value `true`:

```
[calculator.parameters]
xc:                            pw-lda
vdw_correction_hirshfeld:      true
use_gpu:                       true
...
```



### `[calculator.kpoints]` (optional)

#### `density`

Instead of giving a `k_grid` explicitly, FHI-vibes can compute `k_grid` such that the density of kpoints does not fall below this value in $\require{mediawiki-texvc} \AA^{-3}$ . This is optional, including `k_grid` in `[calculator.parameters]` is equally valid.

### `[calculator.basissets]`

Specify which basissets to use.

#### `default`

The default basis set to use, can be `light`, `intermediate`,`tight`, or `really_tight`.

#### `fallback`

The fallback option in case the specified basis set could not be found (`intermediate` basis sets are currently not compiled for each element)

#### Species dependent

The basis set can be given per chemical species by including the species and its desired basis set (uncomment, e.g., `O` in the example above.)

### `[calculator.socketio]` (optional)

Set up socket communication via [`SocketIOCalculator`](https://wiki.fysik.dtu.dk/ase/ase/calculators/socketio/socketio.html?highlight=socketio#ase.calculators.socketio.SocketIOCalculator). This has the potential to speed up calculations since a complete restart of FHI-aims after each completed SCF cycle is avoided. This feature is optional but recommended to use when performing calculations for related structures, e.g., during molecular dynamics simulations or phonon calculations.

#### `port`

The socket port to use.

- `null`: don't use the socket.
- `1024`-`65535`: use this port.

