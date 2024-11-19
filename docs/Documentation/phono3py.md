
vibes supports Boltzmann transport equation (BTE) calculations of thermal transport with 
the force constants calculated from the finite differences method. 
The calculation is enabled by interfacing with [Phono3py](https://phonopy.github.io/phono3py/) packages by setting up a `phono3py.in` file. 
A minimal `phono3py.in` would look like

```fo
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4

[phono3py]
supercell_matrix:              [2, 2, 2]
```

for performing a phono3py calculation for the structure in `geometry.in` with a Lennard-Jones calculator and a $2 \times 2 \times 2$  supercell.

??? info "Click: Default values for the complete list of supported keywords"
    ```
    [phono3py]
    supercell_matrix:              [1, 1, 1]
    displacement:                  0.03
    cutoff_pair_distance:          100.0
    symprec:                       1e-05
    is_diagonal:                   True
    is_plusminus:                  auto
    log_level:                     2
    workdir:                       phono3py
    ```

## Sections

### `[phono3py]`
#### `supercell_matrix`
`list`: A $3 \times 1$ or $9 \times 1$ array specifying the [supercell matrix](../Tutorial/2_phonopy.md#supercell-matrix-supercell_matrix) used for setting up the supercell.

#### `displacement`

`float`: the distance in Å used for the finite displacement.

#### `cutoff_pair_distance`

`float`: cutoff distance for 3rd-order force constants. 
Useful for reducing displaced supercells for 3rd-order force constants.
It corresponds to `phono3py` settings tag [`CUTOFF_PAIR_DISTANCE`](https://phonopy.github.io/phono3py/command-options.html#cutoff-pair-or-cutoff-pair-distance-cutoff-pair-distance)

#### `symprec`

`float`: tolerance for symmetry detection

#### `is_diagonal`

`True/False`: corresponds to the  `phonopy` or `phono3py` settings tag [`DIAG`](https://phonopy.github.io/phonopy/setting-tags.html#diag)

#### `is_plusminus`

`True/False/auto`: corresponds to the  `phonopy` or `phono3py` settings tag [`PM`](https://phonopy.github.io/phonopy/setting-tags.html#pm)

#### `log_level`

`int`: Verbosity control. Can be 0, 1, or 2.

#### `workdir`

`str`: The working directory for running the calculations and storing output files. Default is `phono3py`.

