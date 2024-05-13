!!! info
	A hands-on example for setting up and running a molecular dynamics run [Tutorial](../Tutorial/3_md_intro.md).

vibes supports running molecular dynamics simulations via [ASE](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md) in NVE, NVT, and NPT ensembles. A minimal `md.in` for running 1000 steps of [Velocity Verlet dynamics](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.verlet) with a Lennard Jones calculator would be

```fo
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4

[md]
driver:                        VelocityVerlet
timestep:                      1
maxsteps:                      1000
```

??? info "Click: Default values for the complete list of supported keywords"
	NVE dynamics using Velocity Verlet propagation (`vibes template md`)
    ```
    [md]
    driver:                        VelocityVerlet
    timestep:                      1
    maxsteps:                      1000
    compute_stresses:              False
    workdir:                       md

    [md.kwargs]
    logfile:                       md.log
    ```
    NVT ensemble using a [Langevin](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin) thermostat (`vibes template md --nvt`):
    ```
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
    NPT ensemble using a [Berendsen](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.nptberendsen) thermostat and barostat (`vibes template md --npt`):
    ```
    [md]
    driver:                        NPTBerendsen
    timestep:                      1
    maxsteps:                      1000
    compute_stresses:              False
    workdir:                       md

    [md.kwargs]
    temperature:                   300
    taut:                          500.0
    taup:                          1000.0
    pressure:                      1.01325
    compressibility:               4.57e-05
    logfile:                       md.log
    inhomogeneous:                 False
    ```


## Sections

### `[md]`
#### `driver`
`str`:  The MD algorithm. Either `VelocityVerlet` (NVE ensemble), `Langevin` (NVT ensemble), or `NPTBerendsen` (NPT ensemble).

#### `timestep`

`float`: the integration timestep in fs.

#### `maxsteps`

`int`:  the number of timesteps to perform.

#### `compute_stresses`

`bool/int`: specify whether to compute stress during the MD simulation (more costly in _ab initio_ MD).

#### `workdir`

`str`: the working directory for the MD simulation.

### `[md.kwargs]`
These are keyword arguments  that go straight to the ASE class implementing the MD algorithm, e.g.,

```python
cls = Langevin

md = cls(**settings.md.kwargs)
```

The keywords are documented in ASE:

- [`VelocityVerlet`](https://wiki.fysik.dtu.dk/ase/ase/md.html#ase.md.verlet.VelocityVerlet)
- [`Langevin`](https://wiki.fysik.dtu.dk/ase/ase/md.html#ase.md.langevin.Langevin)
- [`NPTBerendsen`](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.nptberendsen)

For `NPTBerendsen`, we add the `inhomogeneous` keyword, which decides if the `Inhomogeneous_NPTBerendsen` driver is used instead of `NPTBerendsen`. With `Inhomogeneous_NPTBerendsen`, the basis vectors are scaled independently, i.e. the size of the unit cell can change in three directions, but the angles remain constant. By default, `inhomogeneous=False`.

!!! warning
	In ASE, the `temperature` is usually given as an energy in ${\rm eV}$. In FHI-vibes, we use Kelvin consistently. `temperature: 300` thus corresponds to setting the thermostat to $300\,{\rm K}$.
