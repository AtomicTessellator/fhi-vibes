!!! info
	An hands-on example for setting up and running a relaxation can be found in the [Tutorial](../Tutorial/1_geometry_optimization.md).

vibes supports geometry optimization by setting up a `relaxation.in` file. A minimal `relaxation.in` would look like

```fo
[files]
geometry:                      geometry.in

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.4

[relaxation]
driver:                        BFGS
fmax:                          0.001
```

for performing a BFGS optimization of the structure found in `geometry.in` until forces are converged below $\require{mediawiki-texvc} 1\,\text{meV}/\AA$.

??? info "Click: Default values for the complete list of supported keywords"
    ```
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

## Sections

### `[relaxation]`
Instructions to set up a geometry optimization workflow using an [ASE optimizer class](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-ase.optimize).

#### `driver`
Currently only [BFGS](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#bfgs) is supported, which is Quasi-Newton method using the [Broyden–Fletcher–Goldfarb–Shanno algoritm](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) to obtain an estimation of the Hessian.

#### `fmax`

`float`: Maximum residual force in $\require{mediawiki-texvc} \text{eV}/\AA$ (for the stress components: $\text{eV}/\AA^3$ ).

#### `unit_cell`

`True/False`: relax the unit cell using [`ase.constraints.ExpCellFilter`](https://wiki.fysik.dtu.dk/ase/ase/constraints.html?highlight=expcellfilter#ase.constraints.ExpCellFilter)

#### `fix_symmetry`

`True/False`: keep the spacegroup of the system fixed using [`ase.constraints.FixSymmetry`](https://wiki.fysik.dtu.dk/ase/dev/ase/constraints.html?highlight=fixsymmetry#ase.spacegroup.symmetrize.FixSymmetry)

#### `hydrostatic_strain`

`True/False`: apply isotropic pressure instead of stress for cell deformation, see [here](https://wiki.fysik.dtu.dk/ase/ase/constraints.html?highlight=expcellfilter#ase.constraints.ExpCellFilter)

#### `constant_volume`

`True/False`: keep volume constant, see [here](https://wiki.fysik.dtu.dk/ase/ase/constraints.html?highlight=expcellfilter#ase.constraints.ExpCellFilter)

#### `scalar_pressure`

`float`: apply external pressure given in $\text{eV}/\AA^3$

#### `mask`

`list`: `(6, 1)` shaped mask to enable/disable relaxation of strain components in Voigt notation, e.g., to allow for shape-preserving relation. See [ASE documentation of `mask` keyword in `ExpCellFilter` class for details](https://wiki.fysik.dtu.dk/ase/ase/constraints.html?highlight=expcellfilter#ase.constraints.ExpCellFilter). **Using the `mask` keyword will alter the effective stress used to relax the lattice and changes the behavior of `fmax`. Check!**.

#### `decimals`

`int`: number of digits used to round positions before writing `geometry.in.next_step`

#### `symprec`

`float`: symmetry precision used for detecting space group when `fix_symmetry=True`.

#### `workdir`

The working directory for running the calculations and storing output files.

### `[relaxation.kwargs]`

These keywords are used 1:1 by the ASE optimizer class, e.g.

```pytho
cls = ase.optimize.BFGS

optimzer = cls(**settings.relaxation.get("kwargs"))
```

#### `maxstep`

`float`: largest allowed move

#### `logfile`

`str`: logfile for the relaxation

#### `restart`

`str`: use this file to store restart information