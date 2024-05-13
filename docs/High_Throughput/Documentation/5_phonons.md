```
[phonopy]
supercell_matrix = [2, 2, 2]
walltime = 3500
serial = True
displacement = 0.01

[phonopy.convergence]
minimum_similarity_score = 0.80
sc_matrix_base = [1, 1, 1]

[gruneisen]
volume_factors = [0.99, 1.01]

```
## Sections

### `[phonopy]`

Parameters for phonopy calculations. Most keywords are the same. For full documentation of those see the `phonopy` section in the documentation

#### `serial`

`bool`: If True use serial calculations instead of parallel (calculate all supercells in one calculation v. separately) (Default is True)

#### `convergence`

`bool`: If True do phonon supercell convergence with the defaults defined in `phonopy.convergence` section (Default is False)

### `[phonopy.convergence]`

Section used to define the phonon convergence parameters. If both sets of defaults are desired then

#### `minimum_similarity_score`

`float`: Minimum Tanimoto similarity score to consider the phonon calculations converged with respect supercell size (Default is 0.80)

#### `sc_matrix_base`

`list(int)`: Base supercell matrix use to increase supercell size (Default is `phonopy.supercell_matrix`).

In the example above the next supercell matrix tested would be `[3, 3, 3]` without it the next would be `[4, 4, 4]`

### `[gruneisen]`

Set up to calculate the gruneisen parameters from finite difference using phonopy

#### `volume_factors`

`list(float)`: ratio of the equilibrium volume to calculate the phonons for.
