```
[relaxation]
use_ase_relax = False

[relaxation.1]
basis = light
method = trm
fmax = 1e-3
relax_unit_cell = full

[relaxation.2]
basis = intermediate
method = trm
fmax = 1e-3
relax_unit_cell = full

.
.
.

[relaxation.n]
basis = basis_set
method = relaxation_method
fmax = 1e-3
relax_unit_cell = full
```

## Sections

### `[relaxation]`

Sections to do relaxation of structures. This is a general definition for various step

#### `use_ase_relax`

`bool`: True if you want to use the ASE relaxation defined in the `relaxation` portion of the documentation

### `[relaxation.1]`

The first step of the relaxation. If `use_ase_relax` is true replace these keywords with those from the [`relaxation`](../../../Documentation/relaxation) section of the non-high throughput workflows. Otherwise use the electronic structure relaxation keywords. These are examples for `FHI-aims`

#### `basis`

`str`: keyword for the basis set to use for the relaxation step

#### `method`

`str`: Relaxation method used for the calculation

#### `fmax`

`float`: Maximum residual force before ending the relaxation

#### `relax_unit_cell`

`str`: How to relax the unit cell within `FHI-aims` either `full`, `fixed_angles` or `none`

### `[relaxation.2]`

The second step of the relaxation (parameters same as the first step)

### `[relaxation.n]`

The n<sup>th</sup> step of the relaxation

