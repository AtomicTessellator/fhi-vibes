```
[section.stop_if]
condition_list:     ["max(frequencies) < 10.0", ["frequencies", "max", "lt", 10.0, [0.5, 0.0, 0.0]]]
external_functions: ["stop_conditions_module.function_name"]
```

## `[stop_if]`

### `condition_list`:

`list(str or list)`: A list of all conditions to use to defuse the workflow. For string entries the format must be: `"{vector_op}({property}) {comparison_op} {value}"` or `"{property} {comparision_op} {value}"`. For lists the format must be: `["{property}", "{vector_op}", "{comparison_op_name}, {value}, {optional_kwarg_1}, ..., {optional_kwarg_n}]`. For properties with no `optional_kwargs` those elements should not be used.

The available `vector_ops` are an empty string or:

- `max`
- `min`
- `mean`
- `median`
- `std`

The available `comparison_op : comparison_op_names` are:

- `== : eq`
- `!= : ne`
- `<= : le`
- `>= : ge`
- `< : lt`
- `> : gt`

The available properties are defined for each task.

For single-point calculations and relaxations the allowed properties are defined by ASE and use its standards units. The full list of them is defined here:

- `volume`
- `energy`
- `free_energy`
- `lattice_parameters`
- `lattice_angles`

Additionally, for calculations using the `Aims` calculator the estimated `bandgap` (in eV) can be used.

For `phonopy` calculations the allowed properties all use its units and are listed with the following format `{property} ({optional_kwarg}={default_value})`:

- `heat_capacity (T=300)`
- `free_energy (T=300)`
- `entropy (T=300)`
- `frequencies (q=[0, 0, 0])`
- `theta_D`
- `theta_D_infty`
- `theta_P`
- `n_atoms`

The Debye temperature approximations `theta_D`, `theta_D_infty`, and `theta_P` are all defined in `vibes.phonopy.wrapper.get_debye_temperature`. If a string is used for a property with an `optional_kwarg` then the `default_value` will be used.

For `md` and `statistical_sampling` calculations the only allowed property is `sigma`.

### external_functions

`list(str)`: A list of all function paths to user-defined stopping functions. To be used by the workflows the submodule of all functions must be inside the `PYTHONPATH` environment variable

