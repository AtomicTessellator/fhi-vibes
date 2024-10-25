# Thermal lattice expansion

## Recap

**Thermal expansion** is the volume change of materials in response of a change in temperature. For an ideal harmonic system, the harmonic Hamiltonian is independent of the lattice parameters, thus, the lattice expansion vanishes. To correctly estimate lattice expansion, it is essential to consider anharmonic effect. 
One way of estimating the lattice expansion is through Quasi-Harmonic Approximation. You might have already learnt that in [FHI-aims tutorials (Exercise 4: Lattice Expansion in the Quasi-Harmonic Approximation)](https://fhi-aims-club.gitlab.io/tutorials/phonons-with-fhi-vibes/phonons/4_QHA/exercise-4/). Here we are introducing an alternative way of estimate the thermal lattice expansion using NVT molecular dynamics simulations.

## Thermal lattice expansion with MD

In the geometry optimization, the lattice parameter are relaxed at equilibrium (0 K) to remove the internal forces and pressure. However, at finite temperatures, the atoms move around equilibrium, which leads to instantaneous pressure. [As discussed earlier](3_md_intro.md), the pressure decouples into kinetic and potential contribution, $\left\langle p(T) \right\rangle = \left\langle p_{\rm Kin}(T) \right\rangle+ \left\langle p_{\rm Pot}(T) \right\rangle~,$ that can be calculated from the MD simulation. Relaxation of the lattice parameters under the inverse of the expectation values of the pressure, $- \left\langle p(T) \right\rangle$, will remove the residual pressure at finite temperature, $T$. 

## Inspect the pressure using `FHI-vibes`

We use the Lennard-Jones Argon Rapid Prototyping as an example.

After we finish MD at 20 K. We can inspect the pressure from the trajectory file 
```
vibes info md trajectory.nc
```
This command tells you the pressure from the trajectory

```
[info]              Summarize Potential Pressure
Simulation time:                 20.000 ps (5001 of 5001 steps)
Pressure:                          0.052705 +/-   0.005952 GPa
Pressure (last 1/2):               0.053385 +/-   0.005134 GPa
Pressure (last 1/2):               0.000333 +/-   0.000032 eV/AA**3
[info]              Summarize Total Pressure (Kinetic + Potential)
Simulation time:                 20.000 ps (5001 of 5001 steps)
Pressure:                          0.060044 +/-   0.006160 GPa
Pressure (last 1/2):               0.060815 +/-   0.005197 GPa
Pressure (last 1/2):               0.000380 +/-   0.000032 eV/AA**3
```

We take the `Pressure (last 1/2)` as expectation value of the internal pressure at 20 K. And **relax the primitive cell** under the negative of this pressure. We recommend to relax the primitive cell with `fix_symmetry: True` to reserve the space group of the lattice. Here is an example for `relaxation.in`

??? example "Example of `relaxation.in`"
    ```
    [files]
    geometry:                      geometry.in

    [calculator]
    name:                          lj

    [calculator.parameters]
    sigma:                         3.405
    epsilon:                       0.010325
    rc:                            8.0


    [relaxation]
    driver:                        BFGS
    fmax:                          0.0001
    unit_cell:                     True
    fix_symmetry:                  True
    hydrostatic_strain:            False
    constant_volume:               False
    scalar_pressure:               -0.000380
    decimals:                      12
    symprec:                       1e-05
    workdir:                       relaxation

    [relaxation.kwargs]
    maxstep:                       0.2
    logfile:                       relaxation.log
    restart:                       bfgs.restart
    ```

After the relaxation under the negative external scalar pressure. We use the relaxed structure  `geometry.in.next_step` to repeat the MD calculation at 20 K, for checking the residual pressure at this lattice constant.
With the command `vibes info md trajectory.nc`,  we can inspect the pressure of this MD trajectory.

```
[info]              Summarize Total Pressure (Kinetic + Potential)
Simulation time:                 20.000 ps (5001 of 5001 steps)
Pressure:                         -0.041640 +/-   0.005653 GPa
Pressure (last 1/2):              -0.040865 +/-   0.003994 GPa
Pressure (last 1/2):              -0.000255 +/-   0.000025 eV/AA**3
```

The pressure shows  $-0.000255 ~\rm{eV/\mathring{A}^3}$, which indicates that the unit cell is over-expanded after the first iteration. Then, we **subtract this external pressure on top of the first pressure**, $-0.000380 - -0.000255 = -0.000125$ for the scalar pressure for the next iteration. After the second iteration of primitive cell relaxation at `scalar_pressure: -0.000125` and MD simulation at 20 K. The pressure shows:

```
[info]              Summarize Total Pressure (Kinetic + Potential)
Simulation time:                 20.000 ps (5001 of 5001 steps)
Pressure:                         -0.003424 +/-   0.005817 GPa
Pressure (last 1/2):              -0.003018 +/-   0.004259 GPa
Pressure (last 1/2):              -0.000019 +/-   0.000027 eV/AA**3
```

Now the residual pressure $-0.000019 ~\rm{eV/\mathring{A}^3}$ is within the standard error.

## Helper

An example is prepared in [FHI-vibes Gitlab repository](https://gitlab.com/vibes-developers/vibes) at `examples/lattice_expansion/lj_argon_60K`. 

1. `geometry.in` file is the relaxed primitive cell at equilibrium.
2. Navigate to `MD_00` and edit `Makefile`, make sure `T=60` is your target temperature.
3. Run `make all` to make supercell and run MD at your target temperature.
4. Get the `Pressure (last 1/2)` from the printed information
5. Navigate to `relaxation_01` and edit relaxation.in, make sure `scalar_pressure` is the inverse of the pressure from `MD_00`
6. Relax the lattice with command `vibes run relaxation`
7. Navigate to `MD_01` and repeat the workflow from **Step 2**

## TODO:

Calculate thermal lattice expansion at other temperatures, for example 20, 30, 40, 50, 60 K. And plot primitive cell volume as a function of temperatures.
Calculate lattice expansion using quasi-harmonic method, and compare the results.

??? info "thermal lattice expansion"
    ![image](assets/lattice_expansion.png)

