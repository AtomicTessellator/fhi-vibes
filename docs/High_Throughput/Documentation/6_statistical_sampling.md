```
[statistical_sampling]
phonon_file = path/to/phonopy/trajectory.son
supercell_matrix = [-1,1,1,1,-1,1,1,1,-1]
temperatures = [300, 600]
debye_temp_fact = [1.0]
serial = True
n_samples = 1
plus_minus = True
mc_rattle = False
quantum = True
deterministic = True
zacharias = True
gauge_eigenvectors = True
ignore_negative = False
failfast = True
random_seed = 13
propagate = False
```

## Sections

### `[statistical_sampling]`

Used for Monte Carlo sampling of a system for anharmonicity quantification

#### `phonon_file`

`str`: The trajectory file used for generating the thermally displaced structures

#### `supercell_matrix`:

`list(int)`: The supercell matrix for the calculation, if not given use the one from the phonopy calculation. If the supercell matrix is different from the one in the `phonon_file` the phonopy force constants will be remapped onto the new supercell.

#### `temperatures`:

`list(float)`: list of temperatures to calculate the anharmonicity at

#### `debye_temp_fact`:

`list(float)`: list of multipliers to add temperatures that are factors the materials Debye temperature

#### `serial`:

`bool`: If True then do this in serial (Default is True)

#### `n_samples`:

`int`: number of samples to calculate for each temperature (Default is 1)

#### `plus_minus`:

`bool`: Use the deterministic sampling regime from Zacharias, et al (Default is True)

#### `deterministic`:

`bool`: If True populate all phonon modes to k_B T energy (Default is True)

#### `gauge_eigenvectors`:

`bool`: If True use a plus minus gauge for the eigenmodes (Default is True)

#### `mc_rattle`:

`bool`: If True rattle the structures using a Monte Carlo rattling method (Default is False)

#### `quantum`:

`bool`: If True populate phonon modes according to a Bose-Einstein distribution (Default is False)

#### `ignore_negative`:

`bool`: If True ignore all imaginary modes (Default is False)

#### `failfast`:

`bool`: If True If True fail if any imaginary modes are present or acustic modes are not near zero at Gamma (Default is True)

#### `random_seed`:

`int`: The seed for random number generator (Default is Random number)

#### `propagate`:

`bool`: If True propagate the structure forward in time somewhat with ASE (Default is False)
