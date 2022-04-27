# Singlepoint Calculations

!!! info
	We assume you are familiar with `FHI-aims` calculations and are here to learn how to perform them with `FHI-vibes`. The tutorial is however transferable to any other calculator supported by ASE.

In this tutorial, you will learn how to perform singlepoint calculations with `FHI-vibes`.

## Silicon in Equilibrium

As discussed in the introduction, this requires two files. The geometry of the system given by

??? info "geometry.in"

    ```
    lattice_vector 0.0000000000000000 2.72 2.72
    lattice_vector 2.72 0.0000000000000000 2.72
    lattice_vector 2.72 2.72 0.0000000000000000
    atom_frac 0. 0. 0. Si
    atom_frac .25 .25 .25 Si
    ```

and a file desribing the computational taks, which contains the
calculator settings discussed in the introduction. Additionally,
it contains a section `files` pointing to the geometry that shall
be calculated.

??? info "aims.in"

    ```
    [files]
    geometry:                      geometry.in
    
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

You can run the calculation interactively via

```
vibes run singlepoint aims.in | tee log.aims
```

or by submitting it to a queue on a computing cluster.

??? info "Example `submit.sh` for `slurm` queue manager"
    ```
    #!/bin/bash -l

    #SBATCH -J md|vibes
    #SBATCH -o log/md.%j
    #SBATCH -e log/md.%j
    #SBATCH --mail-type=all
    #SBATCH --mail-user=your@mail.com
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=32
    #SBATCH --ntasks-per-core=1
    #SBATCH -t 24:0:00
    
    # Make sure that the correct Python environment is set up, e.g.
    module load miniconda/3/4.5.4
    source activate vibes
    
    vibes run singlepoint aims.in &> log.aims
    ```

The calculation should only take a few seconds and yield the following output:

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
    [calculator]     compute_forces: True
    [calculator]     k_grid: [8, 8, 8]
    [calculator]     sc_accuracy_rho: 1e-06
    [calculator]     relativistic: atomic_zora scalar
    [calculator]     compensate_multipole_errors: False
    [calculator]     output_level: MD_light
    [calculator]     use_pimd_wrapper: ('localhost', 10011)
    [calculator]     aims_command: /u/christia/Codes/vibes_v2/run_aims.sh
    [calculator]     species_dir: /draco/u/christia/Codes/vibes_v2/tutorials/SPC_2/aims/basissets
    [socketio]     Use SocketIO with host localhost and port 10011
    [backup]       /draco/u/christia/Codes/vibes_v2/tutorials/SPC_2/aims/calculations does not exists, nothing to back up.
    [vibes]        Compute structure 1 of 1: working
    [vibes]        Compute structure 1 of 1: finished.
    ```

The calculation will create a working directory called `aims`, the traditional input files for FHI-aims (control.in and geometry.in) as well as the output file
`aims.out` can be found in `aims/calculations/`. Additionally, vibes produces a trajectory file `aims/trajectory.son`, which contains all salient information for
postprocessing and is particularly useful for the adavanced tasks tackled in the next tutorials. The data from this file can be extracted and stored
in an `xarray.Dataset` in `trajectory.nc`, see [the documentation on output files](../Documentation/output_files.md). For this purpose, run

```
vibes output trajectory aims/trajectory.son
```

which yields the output
??? info "Output of `vibes output trajectory`"
    ```
    Extract Trajectory dataset from <Command trajectory>
    [trajectory]   Parse `aims/trajectory.son`
    [son] read file:  aims/trajectory.son
    [son] process:    |||||||||||||||||||||||||||||||||||||  2/2
    [trajectory]   .. time elapsed: 0.000s
    [trajectory]   .. create atoms
    [trajectory]   .. time elapsed: 0.256s
    [trajectory]   .. done in 0.257s
    [trajectory]   Get positions from trajectory
    * Message from file vibes/trajectory/trajectory.py, line 192, function times:
    --> time unit not found in trajectory metadata, use ase.units.fs

    ** Warning from file vibes/trajectory/trajectory.py, line 198, function times:
    --> no time steps found, return time as index
    
    /u/christia/.local/lib/python3.7/site-packages/numpy/core/fromnumeric.py:3373: RuntimeWarning: Mean of empty slice.
      out=out, **kwargs)
    /u/christia/.local/lib/python3.7/site-packages/numpy/core/_methods.py:170: RuntimeWarning: invalid value encountered in double_scalars
      ret = ret.dtype.type(ret / rcount)
    [trajectory]   .. time elapsed: 0.071s
    [trajectory]   Get velocities from trajectory
    [trajectory]   .. time elapsed: 0.001s
    ** Warning from file vibes/trajectory/trajectory.py, line 540, function set_displacements:
    --> SUPERCELL NOT SET, compute w.r.t to reference atoms
    
    [trajectory]   Compute displacements
    [trajectory]   .. time elapsed: 0.002s
    [trajectory]   Get pressure from trajectory
    [trajectory]   .. time elapsed: 0.001s
    Trajectory dataset written to trajectory.nc
    ```


## Multiple Singlepoint Calculations in One Run

`FHI-vibes` offers the possibility to run a set of related calculations in a single run, where "related calculations" means that the input geometries are allowed to differ in their positions and/or lattice.
The stoichometry, the number of atoms, as well as the computational settings must be the same for all systems. For instance, one can run singlepoint calculations for the following two geometries

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

which only differ by in the first fractional coordinate of the first atom, by using the
`geometries` tag in the file `aims.in`


??? info "aims.in"

    ```
    [files]
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

Note, that `FHI-vibes` supports wildcards to read input files and will sort the input files found by this wildcard alphabetically.


Again, you can run the calculation interactively via

```
vibes run singlepoint aims.in | tee log.aims
```

or by submitting to queue as above.

## Submit calculation on a cluster

To efficiently perform _ab initio_ calculations for systems larger than a few atoms, you will need a workstation or access to a computing cluster. To submit a `vibes` simulation to your cluster, follow these steps:

1. [Install `FHI-vibes` on your cluster](../../#installation),
2. set up a calculation as you have done earlier on your laptop,
3. submit the `vibes run` command to the queue.

??? info "Example `submit.sh` for `slurm` queue manager"

    ```
    #!/bin/bash -l
    
    #SBATCH -J md|vibes
    #SBATCH -o log/md.%j
    #SBATCH -e log/md.%j
    #SBATCH --mail-type=all
    #SBATCH --mail-user=your@mail.com
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=32
    #SBATCH --ntasks-per-core=1
    #SBATCH -t 24:0:00
    
    vibes run singlepoint aims.in
    ```

The log file will be written to a `log` folder.

For default tasks, `vibes` can write and submit the `submit.sh` file in a single command, see [documentation of slurm submission](../Documentation/input_files_slurm.md)

## Restart a calculation

If your calculation does not fit into a walltime or stops for another reason before the total number of simulation steps is reached, you can simply resubmit `vibes run`. It will restart the calculation from the last completed step. This also holds for `vibes run phonopy` and `vibes run md` as explained later.

### Automatic restarts

Optionally, `vibes` can restart the job by itself using a `[restart]` section in the input file, e.g., `aims.in`. To this end, add

```
[restart]
command = sbatch submit.sh
```

to your `aims.in`, where `sbatch submit.sh` is the command you use to submit the calculation to the queue. `vibes` will run this command shortly before the walltime is over to restart the job.

## Energy Unit Warning When Using `socketio`

If you are using a `socketio` calculator for single point calculations then the conversion from Hartree to eV will be done using ASE.
This can lead to slight deviations between the parameters listed in the `trajectory.son` and the output file from the calculation if the base calculator does not use the same conversion factor between Hartree and eV.
The `SocketIOCalculator` will use the CODATA 2014 standard as of  `ASE` version `3.23`.
