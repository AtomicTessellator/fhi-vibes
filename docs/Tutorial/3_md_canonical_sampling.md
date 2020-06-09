# Molecular Dynamics: Canonical Sampling

The aim of this tutorial is to learn how to perform a molecular dynamics simulation in the canonical ensemble and thermalize a system to a target temperature.

!!! warning
	This tutorial mimics the essential steps for performing MD simulations in bulk systems. How you obtain initial structures in your project is, of course, highly dependent on the systems you aim to study etc.

!!! info
	For this tutorial, we study solid Lennard-Jones Argon at $20\,{\rm K}$. All steps are transferable to performing _ab initio_ molecular dynamics simulations by exchanging the calculator with `FHI-aims`.

##  Structure preparation

### Generate a structure

We first create an input structure for solid Argon. To this end, we first generate a primitive cell of cubic Argon, which can be performed with ASE:

```python
>>> from ase.build import bulk
>>> 
>>> atoms = bulk('Ar')
>>> 
>>> atoms.write('geometry.in.primitive', format='aims')
```

From this primitive cell, we can use the CLI to generate a supercell of about 100 atoms:

```
vibes utils make_supercell geometry.in.primitive -n 100
```

this will try to find a cubic-as-possible supercell with roughly 100 atoms and write it to an input file.

??? info "Output  of `vibes utils make_supercell geometry.in.primitive -n 100`"
    ```
    Find supercell for
    [vibes]        Geometry info
      input geometry:    Ar
      Symmetry prec.:    1e-05
      Number of atoms:   1
      Species:           Ar (1)
      Periodicity:       [ True  True  True]

      Spacegroup:          Fm-3m (225)
      Wyckoff positions:   1*a
      Equivalent atoms:    1*0
    
    Settings:
      Target number of atoms: 100
    
    Supercell matrix:
     python:  [-3,  3,  3,  3, -3,  3,  3,  3, -3]
     cmdline: -3 3 3 3 -3 3 3 3 -3
     2d:
    [[-3, 3, 3],
     [3, -3, 3],
     [3, 3, -3]]
    
    Superlattice:
    [[15.78  0.    0.  ]
     [ 0.   15.78  0.  ]
     [ 0.    0.   15.78]]
    
    Number of atoms:  108
      Cubicness:         1.000 (1.000)
      Largest Cutoff:    7.890 AA
    
    Supercell written to geometry.in.primitive.supercell_108
    ```

In this case it should find a perfectly cubic supercell with 108 atoms and write it to `geometry.in.primitive.supercell_108`. You can obtain detailed information about the structure by running

```
vibes info geometry geometry.in.primitive.supercell_108
```

??? info "Output of `vibes info geometry geometry.in.primitive.supercell_108`"
    ```
    [vibes]        Geometry info
      input geometry:    Ar
      Symmetry prec.:    1e-05
      Number of atoms:   108
      Species:           Ar (108)
      Periodicity:       [ True  True  True]
      Lattice:  
        [15.78  0.    0.  ]
        [ 0.   15.78  0.  ]
        [ 0.    0.   15.78]
      Cubicness:         1.000 (1.000)
      Largest Cutoff:    7.890 AA

      Spacegroup:          Fm-3m (225)
      Wyckoff positions:   108*a
      Equivalent atoms:    108*0
    
    Cell lengths and angles [Å, °]:
      a, b, c:     15.7800     15.7800     15.7800
      α, β, γ:     90.0000     90.0000     90.0000
      Volume:             3929.3526 Å**3
      Volume per atom:      36.3829 Å**3
    ```

Additionally you can inspect the generated structure with the structure viewer of your choice, e.g., with [`jmol`](http://jmol.sourceforge.net/).

??? info "`jmol geometry.in.primitive.supercell_108`"
	![image](../assets/geometry.in.supercell.png)



Assuming that we are happy with this structure, we save it as our supercell:

```
mv geometry.in.primitive.supercell_108 geometry.in.supercell
```

### Pre-thermalize the structure

To speed up the thermalization, we can pre-thermalize the system by giving momenta to the atoms according to a Maxwell-Boltzmann distribution at our target temperature of $ 20\,{\rm K}$. This can be done with the CLI utility `create-samples`:

```
vibes utils create-samples geometry.in.supercell -T 20
```

??? info "Output of `vibes utils create-samples geometry.in.supercell -T 20`"
    ```
    vibes CLI: create_samples
    [vibes]        Geometry info
      input geometry:    Ar
      Symmetry prec.:    1e-05
      Number of atoms:   108
      Species:           Ar (108)
      Periodicity:       [ True  True  True]

      Spacegroup:          Fm-3m (225)
      Wyckoff positions:   108*a
      Equivalent atoms:    108*0
    [vibes]        Geometry info
      input geometry:    Ar
      Symmetry prec.:    1e-05
      Number of atoms:   108
      Species:           Ar (108)
      Periodicity:       [ True  True  True]
    
      Spacegroup:          Fm-3m (225)
      Wyckoff positions:   108*a
      Equivalent atoms:    108*0
    [vibes]        Use Maxwell Boltzamnn to set up samples
    [vibes]        Sample   0:
    [vibes]        .. temperature before cleaning:    18.119K
    [vibes]        .. remove net momentum from sample and force temperature
    [vibes]        .. temperature in sample 0:        20.000K
    [vibes]        Sample   0:
    [vibes]        .. temperature in sample 0:        20.000K
    [vibes]        .. written to geometry.in.supercell.0020K
    ```
The geometry written to `geometry.in.supercell.0020K` will now include the appropriate velocities.

We will use this structure and the chosen velocities as the initial structure for the MD run. We suggest to rename this file to `geometry.in` accordingly.

## Run MD

### Prepare `md.in`

Before we can run the MD, we need to create an input file. We can use the CLI command `template` to do this:

```
vibes template lj > md.in 
vibes template md --nvt >> md.in
```

??? info "The generated `md.in`"

    ```
    [files]
    geometry:                      geometry.in
    
    [calculator]
    name:                          lj
    
    [calculator.parameters]
    # parameters for LJ Argon
    sigma:    3.405
    epsilon:  0.010325 
    rc:       13.0


    [md]
    driver:                        Langevin
    timestep:                      1
    maxsteps:                      1000
    compute_stresses:              False
    workdir:                       md
    
    [md.kwargs]
    temperature:                   300
    friction:                      0.02
    logfile:                       md.log
    ```

We suggest to add and/or adjust the following  keywords:

```
[files]
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell

[calculator.parameters]
rc:       8.0

[md]
timestep:                      4
maxsteps:                      7500

[md.kwargs]
temperature:                   20
```

Decreasing `rc` will speed up the calculation, as this is a tutorial. `timestep` can be increased to $4\,{\rm fs}$ for Argon at $20\,{\rm K}$. With `maxsteps: 7500` we will run a total of 7500 MD steps, i.e., $30\,{\rm ps}$ simulation time. **The total simulation time depends on the system and the quantitiy of interest!**`temperature` should be set to $20\,{\rm K}$, our target temperature.
Adding `primitive: geometry.in.primitive` and `supercell: geometry.in.supercell` in the `[files]` section is not necessary to run the calculation. However, `vibes` will automatically attach this information to the trajectory so that it cannot get lost. This also makes life easer when post processing. For example, the displacements $\Delta {\bf R}_I(t)$ can only be properly calculated, when the reference supercell is known.

The final `md.in` should look like this:

```
[files]
geometry:                      geometry.in
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell

[calculator]
name:                          lj

[calculator.parameters]
# parameters for LJ Argon
sigma:    3.405
epsilon:  0.010325 
rc:       8.0


[md]
driver:                        Langevin
timestep:                      4
maxsteps:                      7500
compute_stresses:              False
workdir:                       md

[md.kwargs]
temperature:                   20
friction:                      0.02
logfile:                       md.log
```

We are now ready to  run the simulation!

### Run the calculation

You can run this calculation with the CLI command `run`. We recommend to save it's output, e.g., with `tee`:

```
vibes run md | tee md.log
```

Depending on you computer, the calculation will take a few minutes.

## Create trajectory dataset and inspect the simulation

### Process the calculation

The data obtained at each time step will be written to the trajectory file `md/trajectory.son`. The CLI provides a tool to process the trajectory and create an `xarray.Dataset` from it. To this end, run

 ```
vibes output md md/trajectory.son
 ```

This command will create `trajectory.nc`, a dataset representation of the data contained in the MD trajectory saved as an `xarray.Dataset` to a NetCDF file. The included data can be viewed with

```
vibes info netcdf trajectory.nc
```

??? info "Output of `vibes info netcdf trajectory.nc`"

    ```
    <xarray.Dataset>
    Dimensions:              (I: 108, a: 3, b: 3, time: 2501)
    Coordinates:
      * time                 (time) float64 0.0 4.0 8.0 ... 9.996e+03 1e+04
    Dimensions without coordinates: I, a, b
    Data variables:
        positions            (time, I, a) float64 ...
        displacements        (time, I, a) float64 ...
        velocities           (time, I, a) float64 ...
        momenta              (time, I, a) float64 ...
        forces               (time, I, a) float64 ...
        energy_kinetic       (time) float64 ...
        energy_potential     (time) float64 ...
        stress               (time, a, b) float64 ...
        stress_kinetic       (time, a, b) float64 ...
        stress_potential     (time, a, b) float64 ...
        temperature          (time) float64 ...
        cell                 (time, a, b) float64 ...
        positions_reference  (I, a) float64 ...
        lattice_reference    (a, b) float64 ...
        pressure             (time) float64 ...
        pressure_kinetic     (time) float64 ...
        pressure_potential   (time) float64 ...
    Attributes:
        name:             trajectory
        system_name:      Ar
        natoms:           108
        time_unit:        fs
        timestep:         4.000000000000006
        nsteps:           2500
        symbols:          ['Ar', 'Ar', 'Ar', 'Ar', 'Ar', 'Ar', 'Ar', 'Ar', 'Ar', ...
        masses:           [39.948 39.948 39.948 39.948 39.948 39.948 39.948 39.94...
        atoms_reference:  {"pbc": [true, true, true],\n"cell": \n[[ 1.57800000000...
        volume:           3929.352552000002
        raw_metadata:     {"MD": {\n  "type": "molecular-dynamics",\n  "md-type":...
        hash:             0eff05aa63cd4019927c42af74bb0ff0a0e21009
    ```

### View simulation statistics

To get information about the simulation, you can use the CLI command `info md`, which summarizes the simulation and can produce an overview plot:

```
vibes info md trajectory.nc -p
```

This command should tell you, among other things, that the temperature is indeed thermalized to approximately $ 20\,{\rm K}$:

```
...
[info]         Summarize Temperature
Simulation time:            30.000 ps (7501 steps)
Temperature:                    19.826 +/-       1.7902 K
Temperature (1st 1/3):          19.402 +/-       2.2027 K
Temperature (2st 1/3):          20.426 +/-       1.4744 K
Temperature (3st 1/3):          19.651 +/-       1.4218 K
Temperature (last 1/2):         19.863 +/-       1.4301 K
...
```

The pdf file `md_summary.pdf` provides a visualization of the simulated properties for quick sanity checking that  the simulation went according to plan:

??? info "`md_summary.pdf`"
	![image](../assets/md_summary.png)
	
### Visualize trajectory
The trajectory can be exported to an `xyz` file for visualizing the atomic motion, e.g., with [`VMD`](https://www.ks.uiuc.edu/Research/vmd/). To this end, run

```
vibes utils trajectory 2xyz trajectory.nc
vmd trajectory.xyz
```

??? info "`vmd trajectory.xyz`"
	![image](../assets/LJ-Argon.gif)
