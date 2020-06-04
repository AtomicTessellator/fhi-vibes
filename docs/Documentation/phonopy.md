!!! info
	An hands-on example for setting up and running a `phonopy` calculation can be found in the [Tutorial](../Tutorial/2_phonopy.md).

vibes supports [phonon calculations with the finite differences method](../Tutorial/2_phonopy_intro.md#phonons-harmonic-vibrations-in-solids) by setting up a `phonopy.in` file. A minimal `phonopy.in` would look like

```fo
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4

[phonopy]
supercell_matrix:              [2, 2, 2]
```

for performing a phonopy calculation for the structure in `geometry.in` with a Lennard-Jones calculator and a $2 \times 2 \times 2$  supercell.

??? info "Click: Default values for the complete list of supported keywords"
    ```
    [phonopy]
    supercell_matrix:              [1, 1, 1]
    displacement:                  0.01
    is_diagonal:                   False
    is_plusminus:                  auto
    symprec:                       1e-05
    q_mesh:                        [45, 45, 45]
    workdir:                       phonopy
    ```

## Sections

### `[phonopy]`
#### `supercell_matrix`
`list`: A $3 \times 1$ or $9 \times 1$ array specifying the [supercell matrix](../Tutorial/2_phonopy.md#supercell-matrix-supercell_matrix) used for setting up the supercell.

#### `displacement`

`float`: the distance in $\require{mediawiki-texvc} \AA$ used for the finite displacement.

#### `is_diagonal`

`True/False`: corresponds to the  `phonopy` settings tag [`DIAG`](https://phonopy.github.io/phonopy/setting-tags.html#diag)

#### `is_plusminus`

`True/False/auto`: corresponds to the  `phonopy` settings tag [`PM`](https://phonopy.github.io/phonopy/setting-tags.html#pm)

#### `symprec`

`float`: tolerance for symmetry detection

#### `q_mesh`

`list`: the q-points mesh used for postprocessing, e.g., density of states.

#### `workdir`

The working directory for running the calculations and storing output files.
