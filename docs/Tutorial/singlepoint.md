# Singlepoint Calculations

!!! info
	We assume you are familiar with `FHI-aims` calculations and are here to learn how to perform them with `FHI-vibes`. The tutorial is however transferable to any other calculator supported by ASE.

In this tutorial, you will learn how to perform singlepoint calculations with `FHI-vibes`. `vibes` offers the possibility to run a set of related calculations in a single run, where "related calculations" means that the input geometries are allowed to differ by changed positions and/or lattice. The calculations will be performed with similiar computational settings.

## Define Inputs

Take the following 4 input files and copy them into a working directory:

??? info "geometry.in"

    ```
    lattice_vector 0.0000000000000000 2.72 2.72
    lattice_vector 2.72 0.0000000000000000 2.72
    lattice_vector 2.72 2.72 0.0000000000000000 
    atom_frac 0. 0. 0. Si
    atom_frac .25 .25 .25 Si
    ```
    
??? info "geometry.in.000"

    ```
    lattice_vector 0.0000000000000000 2.72 2.72
    lattice_vector 2.72 0.0000000000000000 2.72
    lattice_vector 2.72 2.72 0.0000000000000000 
    atom_frac 0.01 0. 0. Si
    atom_frac .25 .25 .25 Si
    ```

??? info "geometry.in.001"

    ```
    lattice_vector 0.0000000000000000 2.72 2.72
    lattice_vector 2.72 0.0000000000000000 2.72
    lattice_vector 2.72 2.72 0.0000000000000000 
    atom_frac 0.02 0. 0. Si
    atom_frac .25 .25 .25 Si
    ```

??? info "aims.in"

    ```
    [files]
    geometry:                      geometry.in
    geometries:                    geometry.in.???

    [calculator]
    name:                          aims
    socketio:                      true

    [calculator.parameters]
    xc:                            pw-lda
    compute_forces:                true

    [calculator.kpoints]
    density:                       3

    [calculator.basissets]
    Si:                            light
    ```
    
As you see, the three input structures differ by having the first atom displaced along the first lattice vector in `geometry.in.000` and `geometry.in.001`. The three input files are specified in the input file like this:

```
[files]
geometry:                      geometry.in
geometries:                    geometry.in.???
```

Please note, that `vibes` supports wildcards to read input files and will sort the input files found by this wildcard alphabetically. The tag `geometry` specifies the first geometry to be read, which will also serve as a reference geometry. `geometries` specifies all other files with structures that will be calculated afterwards. Note that `geometry` can be omitted if there is no reference structure to be computed.

The calculator is set up from the [`[calculator]` sections in the input file](../Documentation/calculator_setup.md):

```
[calculator]
name:                          aims
socketio:                      true

[calculator.parameters]
xc:                            pw-lda
compute_forces:                true

[calculator.kpoints]
density:                       3

[calculator.basissets]
Si:                            light
```

Let's walk through the settings once:

- `[calculator]`
	- `name: aims` means that `FHI-aims` will be used as explained [here](../Documentation/calculator_setup.md#calculator).
	- `socketio: true` means that the socketIO communication will be used. This will speed up the computation of related structures.
- `[calculator.parameters]`: these are settings that go directly to `control.in`
	- `xc: pw-lda` means that the pw-LDA exchange-correlation functional will be used.
	- `compute_forces: true` means that forces will be computed.
- `[calculator.kpoints]`: this is an optional way of setting k-point grids based on a target density (specifying `k_grid: X Y Z` in `[calculator.parameters]` is also possible!)
	- `density: 3` use a k-point density of at least 3 per $\require{mediawiki-texvc} \AA^{-3}$.
- `[calculator.basissets]`: Details on which basissets to use
	- `Si: light`: use _light default_ basis sets for silicon.

We are good to go!

## Run the calculation

You can run the calculation via

```
vibes run singlepoint aims.in | tee log.aims
```

The calculation should only take a few seconds.

??? info "`log.aims`"

    ```
    [vibes.run]    run singlepoint calculations with settings from aims.in

    [calculator]   Update aims k_grid with kpt density of 3 to [8, 8, 8]
    [calculator]   .. add `sc_accuracy_rho: 1e-06` to parameters (default)
    [calculator]   .. add `relativistic: atomic_zora scalar` to parameters (default)
    [calculator]   .. add `compensate_multipole_errors: False` to parameters (default)
    [calculator]   .. add `output_level: MD_light` to parameters (default)
    [calculator]   Add basisset `light` for atom `Si` to basissets folder.
    [calculator]   Calculator: aims
    [calculator]   settings:
    [calculator]     xc: pw-lda
    [calculator]     k_grid: [8, 8, 8]
    [calculator]     sc_accuracy_rho: 1e-06
    [calculator]     relativistic: atomic_zora scalar
    [calculator]     compensate_multipole_errors: False
    [calculator]     output_level: MD_light
    [calculator]     use_pimd_wrapper: ('localhost', 10001)
    [calculator]     aims_command: run_aims
    [calculator]     species_dir: /home/knoop/local/hilde/examples/tutorial/0_singlepoint_calculations/aims/basissets
    [watchdog]     seems we are not on a cluster, nothing to do for watchdog
    [socketio]     Use SocketIO with host localhost and port 10001
    [backup]       /home/knoop/local/hilde/examples/tutorial/0_singlepoint_calculations/aims/calculations does not exists, nothing to back up.
    [vibes]        Compute structure 1 of 3: working
    [vibes]        Compute structure 1 of 3: finished.
    [vibes]        Compute structure 2 of 3: working
    [vibes]        Compute structure 2 of 3: finished.
    [vibes]        Compute structure 3 of 3: working
    [vibes]        Compute structure 3 of 3: finished.
    [vibes]        done.
    ```
    
The calculation will create a working directory called `aims` and produce a trajectory file `aims/trajectory.son`. The data from this file can be extracted with

```
vibes output trajectory aims/trajectory.son
```

which will create an `xarray.Dataset` in `trajectory.nc`, see [the documentation on output files](../Documentation/output_files.md). The file contains total energies and forces. Additionally, you find the traditional `aims.out` file in `aims/calculations/aims.out`.