<a name="1_GeometryOptimization"></a>

In this tutorial, you will learn how to perform a geometry optimization with `FHI-vibes`.

!!! info
	We give explicit references for LDA-Silicon. When using LJ-Argon, the only difference lies the definition of the calculator in the `[calculator]` section, and the respective structure defined in `geometry.in`.

## Define Inputs

For starting the relaxation, we use the `geometry.in` file for Silicon in the primitive unit cell  discussed in the [introduction](0_intro.md#test-systems) and copy it to the actual (empty) directory in which we are working. Generate a task input file for running a relaxation by copying the [calculator information for your test system](0_intro.md#test-systems) to a file called `relaxation.in`. Next, use the command line interface (CLI) of `FHI-vibes` to obtain default settings for performing the relaxation and appending them to the input file:

```
vibes template relaxation >> relaxation.in
```

In case of LDA-Silicon with `FHI-aims` calculator, the newly generated input file `relaxation.in` should look like this:

??? info "`relaxation.in`"
    ```
	[files]
	geometry:                      geometry.in
	
	[calculator]
    name:                          aims
    socketio:                      True
    
    [calculator.parameters]
    xc:                            pw-lda
    
    [calculator.kpoints]
    density:                       2
    
    [calculator.basissets]
    default:                       light
    
    [relaxation]
    driver:                        BFGS
    fmax:                          0.001
    unit_cell:                     True
    fix_symmetry:                  False
    hydrostatic_strain:            False
    constant_volume:               False
    scalar_pressure:               0.0
    decimals:                      12
    symprec:                       1e-05
    workdir:                       relaxation
    
    [relaxation.kwargs]
    maxstep:                       0.2
    logfile:                       relaxation.log
    restart:                       bfgs.restart
    ```

The settings file template you just generated contains all the necessary settings to set up and run a geometry optimization with `FHI-vibes` using `FHI-aims` as the force/stress calculator. 
`FHI-vibes` will perform a [BFGS optimization of the structure as implemented in ASE](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#bfgs).
The keywords are explained in the [documentation](../Documentation/relaxation.md).

## Run calculation
You can start an interactive calculation with `vibes run relaxation` or by incorporating this command
in the respectives submission file (see Sec. Singlepoint).
We suggest pipe the output, e.g., like this:

```
vibes run relaxation > log.relaxation 
```

`FHI-vibes` will create a working directory with the default name `relaxation` and will handle running the `FHI-aims` calculations necessary to perform the geometry optimization.
The log should read like that: 

??? info "`relaxation.log`"
    ```
      [vibes.run]    run relaxation workflow with settings from relaxation.in

      [relaxation]   ** /draco/u/christia/Codes/vibes_v2/tutorials/GR/relaxation/trajectory.son does not exist, nothing to prepare
      [calculator]   Update aims k_grid with kpt density of 3 to [8, 8, 8]
      [calculator]   .. add `sc_accuracy_rho: 1e-06` to parameters (default)
      [calculator]   .. add `relativistic: atomic_zora scalar` to parameters (default)
      [calculator]   .. add `compensate_multipole_errors: False` to parameters (default)
      [calculator]   .. add `output_level: MD_light` to parameters (default)
      [calculator]   Add basisset `light` for atom `Si` to basissets folder.
      [calculator]   Calculator: aims
      [calculator]   settings:
      [calculator]     xc: pw-lda
      [calculator]     compute_forces: True
      [calculator]     k_grid: [8, 8, 8]
      [calculator]     sc_accuracy_rho: 1e-06
      [calculator]     relativistic: atomic_zora scalar
      [calculator]     compensate_multipole_errors: False
      [calculator]     output_level: MD_light
      [calculator]     compute_analytical_stress: True
      [calculator]     use_pimd_wrapper: ('localhost', 10011)
      [calculator]     aims_command: /u/christia/Codes/vibes_v2/run_aims.sh
      [calculator]     species_dir: /draco/u/christia/Codes/vibes_v2/tutorials/GR/relaxation/basissets
      [relaxation]   filter settings:
      [relaxation]     hydrostatic_strain: False
      [relaxation]     constant_volume: False
      [relaxation]     scalar_pressure: 0.0
      [relaxation]   driver: BFGS
      [relaxation]   settings:
      [relaxation]     type: optimization
      [relaxation]     optimizer: BFGS
      [relaxation]     maxstep: 0.2
      [socketio]     Use SocketIO with host localhost and port 10011
      [relaxation]   filter settings:
      [relaxation]     hydrostatic_strain: False
      [relaxation]     constant_volume: False
      [relaxation]     scalar_pressure: 0.0
      [relaxation]   Start step 0
      [relaxation]   Step 0 finished.
      [relaxation]   .. residual force:  0.000 meV/AA
      [relaxation]   .. residual stress: 289.641 meV/AA**3
      [vibes]        .. Space group:     Fd-3m (227)
      [relaxation]   clean atoms before logging
      [relaxation]   .. log
      [relaxation]   Step 1 finished.
      [relaxation]   .. residual force:  0.000 meV/AA
      [relaxation]   .. residual stress: 3.463 meV/AA**3
      [vibes]        .. Space group:     Fd-3m (227)
      [relaxation]   clean atoms before logging
      [relaxation]   .. log
      [relaxation]   Step 2 finished.
      [relaxation]   .. residual force:  0.000 meV/AA
      [relaxation]   .. residual stress: 0.039 meV/AA**3
      [vibes]        .. Space group:     Fd-3m (227)
      [relaxation]   clean atoms before logging
      [relaxation]   .. log
      [relaxation]   Relaxation converged.
      [relaxation]   done.
    ```

You will find the FHI-aims in- and output in `relaxation/calculation/`, the final converged structure in `relaxation/geometry.in.next_step`, and a summary of the relaxtion path in `relaxation/relaxation.log`.

For a detailed summary of the relaxation path, you may run

```
vibes info relaxation relaxation/trajectory.son
```

??? info "Output"
    ```
    Relaxation info for relaxation/trajectory.son:
    fmax:             1.000e+00 meV/AA
    # Step |   Free energy   |   F-F(1)   | max. force |  max. stress |  Volume  |  Spacegroup  |
    #      |       [eV]      |    [meV]   |  [meV/AA]  |  [meV/AA^3]  |  [AA^3]  |              |

        1    -15748.20070140     -0.605222       0.0000         0.3084     39.773   Fd-3m (227)
    --> converged.
    ```
