```
[md]
phonon_file = path/to/phonopy/trajectory.son
temperatures = [300, 600]
supercell_matrix = [1, 1, 1]
{The rest of the MD parameters described in [md]}
```

## Sections

### `[md]`

#### `phonon_file`

`str`: The trajectory file used for generating the thermally displaced structures

#### `supercell_matrix`:

`list(int)`: The supercell matrix for the calculation, if not given use the one from the phonopy calculation. If the supercell matrix is different from the one in the `phonon_file` the phonopy force constants will be remapped onto the new supercell.

#### `temperatures`:

`list(float)`: list of temperatures to calculate the anharmonicity at
