!!! danger
	In this Tutorial, we assume you are already familiar with running *FHI-aims* calculations, and that you have [installed](../README.md#installation) and [configured](../README.md#configuration) `FHI-vibes` successfully.
	
!!! info
	For vibrational studies, it is crucial to use structures that are accurately  relaxed. Before starting with actual phonon calculations, we thus learn how to perform a geometry optimization with *FHI-vibes*

## <a name="1_GeometryOptimization"></a> Optimize Your Geometry

### Define Inputs

As input, we use an fcc-diamond Silicon crystal:
```
# save this as geometry.in
lattice_vector 0.00 2.72 2.72
lattice_vector 2.72 0.00 2.72
lattice_vector 2.72 2.72 0.00

atom_frac 0.00 0.00 0.00 Si
atom_frac 0.25 0.25 0.25 Si
```

??? info "Use `ase.build` to create geometry input files for common materials"
	You can use _ase_ to create a template Silicon crystal input file with just three lines of code:

    ```python
    """run this in a script or in a python shell"""
    
    from ase.build import bulk

    si = bulk("Si")

    si.write("geometry.in", scaled=True)
    ```

Next, we generate an input file for running a relaxation via the command line interface (CLI) of `FHI-vibes`. Since we want to use `FHI-aims` for performing the DFT calculation of forces and stress, we obtain template settings for setting up an `aims` calculator:

```
vibes template aims > relaxation.in
```

Next we add template settings for performing the relaxation:

```
vibes template relaxation >> relaxation.in
```

The newly generated input file `relaxation.in` should look like this:

??? info "`relaxation.in`"
    ```
    [files]
    geometry:                      geometry.in

    [calculator]
    name:                          aims

    [calculator.parameters]
    xc:                            pw-lda
    k_grid:                        [2, 2, 2]

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

        1    -15745.65358368    -50.423536       0.0000        29.8125     43.172   Fd-3m (227)
        2    -15745.65367432    -50.514166       0.0000         2.9539     43.118   Fd-3m (227)
        3    -15745.65368216    -50.522012       0.0000         0.0130     43.112   Fd-3m (227)
    --> converged.
    ```