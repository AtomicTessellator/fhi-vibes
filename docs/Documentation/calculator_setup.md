# Calculator Setup

FHI-vibes can set up any [ASE calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#module-ase.calculators) for performing a calculation by providing the  calculator class `name` and the respective `parameters` in the input file. If a `module` is specified, `vibes` will attempt to import the calculator from that module instead of `ase`. This can be used to work with custom calculators that are not (yet) included in `ase`.


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

For a non-`ase` calculator, this would be:

```
...
[calculator]
name:                          MyCalculator
module:                        mymodule

[calculator.parameters]
a:                             1.23
...
```

`vibes` will then attempt to import `MyCalculator` from `mymodule` and instantiate it with `a=1.23`. 

## Sections

### `[calculator]`

This section specifies which `ase.Calculator` should be set up and how.

#### `name`

The name of the [ASE calculator class name](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#supported-calculators).

Note that for non-`ase` calculators, `name` must be spelled identically to the class name in the module, i.e. typically `CamelCase`.

#### `module` (optional)

If specified, `vibes` will run `from module import name` to obtain the calculator class, instead of importing it from `ase`. 

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

#### `host`

The IP address to access the socket. Default is `localhost` and will only have to be modified for certain architectures (this will likely be made clear in the system documentation).

#### `port`

The socket port to use.

- `null`: don't use the socket.
- `auto`: Automatically select a port that is not currently in use by the `host` or registered in `/etc/services`
- `1024`-`65535`: If available use this port (if it's not already being used). We recommend using `auto` for all calculations.

#### `unixsocket`

Filename for the unix socket. If this is active TCP/IP socket will not be used (not recommended, but maybe necessary on some systems)