<a name="1_GeometryOptimization"></a>

In this tutorial, you will learn how to perform a geometry optimization with `FHI-vibes`.

!!! info
	We give explicit references for LDA-Silicon. When using LJ-Argon, the only difference lies in how the `calculator` is defined and your `geometry.in` input file.

## Define Inputs

[Choose a test system](0_intro.md#test-systems) and copy the geometry information into a file called `geometry.in`. Next, we generate an input file for running a relaxation via the command line interface (CLI) of `FHI-vibes`. To this end, please copy the [calculator information for your test system](0_intro.md#test-systems) to a file called `relaxation.in`.

Next we add template settings for performing the relaxation:

```
vibes template relaxation >> relaxation.in
```

In case of LDA-Silicon with `FHI-aims` calculator, the newly generated input file `relaxation.in` should look like this:

??? info "`relaxation.in`"
    ```
	[calculator]
    name:                          aims
    
    [calculator.parameters]
    xc:                            pw-lda
    
    [calculator.kpoints]
    density:                       2
    
    [calculator.basissets]
    default:                       light
    
    [calculator.socketio]
    port:                          12345
    
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

The settings file template you just generated contains all the necessary settings to set up and run a geometry optimization with `FHI-vibes` using `FHI-aims` as the force/stress calculator (or Lennard-Jones if you're working with Argon).

## Run calculation
You can start the calculation with `vibes run relaxation`. We suggest pipe the output, e.g., like this:

```
vibes run relaxation > log.relaxation &
```

`vibes` will create a working directory with the default name `relaxation` and will handle running the `aims` calculations and using a [straightforward BFGS algorithm implemented in ASE](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#bfgs). You will find the converged structure in `relaxation/geometry.in.next_step`, and a summary of the relaxtion path in `relaxation/relaxation.log`.

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