```
[relaxation]
use_aims_relax = True
method = trm
fmax = 1e-3
relax_unit_cell = full

[relaxation.1]
basis = light


[relaxation.2]
basis = intermediate

.
.
.

[relaxation.n]
basis = really_tight
```

## Sections

### `[relaxation]`

Sections to do relaxation of structures. This is a general definition for various step

#### `use_aims_relax`

`bool`: True if you want to use the aims relaxation instead of the one defined in the [`relaxation`](../../../Documentation/relaxation) portion of the documentation. For this documentation we will only show keywords if this is true. If this is false then use the keywords you'd normally use for a relaxation

#### `basis`

`str`: keyword for the basis set to use for the relaxation step

#### `method`

`str`: Relaxation method used for the calculation

#### `fmax`

`float`: Maximum residual force before ending the relaxation

#### `relax_unit_cell`

`str`: How to relax the unit cell within `FHI-aims` either `full`, `fixed_angles` or `none`

### `[relaxation.1]`

The first step of the relaxation. Only define parameters that are different from default parameters


### `[relaxation.2]`

The second step of the relaxation (parameters same as the first step)

### `[relaxation.n]`

The n<sup>th</sup> step of the relaxation

