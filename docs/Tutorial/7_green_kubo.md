# Green-Kubo workflow

## Introduction

Green-Kubo theory is based on linear response theory, it calculates the thermal conductivity from equilibrium molecular dynamics. Green-Kubo theory relates the thermal conductivity from the integration of the real time autocorrelation function of heat flux $\boldsymbol{J}$,

$$
    \boldsymbol{\kappa}^{\alpha \beta}(T) = \frac{V}{k_B T^2} \int_{0}^{t_0} \langle \boldsymbol{J}^{\alpha}(t) \cdot \boldsymbol{J}^{\beta}(0) \rangle_T \ dt
$$

where $\boldsymbol\kappa^{\alpha \beta}$ is the thermal conductivity tensor element in Cartesian coordinates of $\alpha, \beta$, with the volume of the supercell volume $V$, Boltzmann constant $k_B$ and temperature $T$.
The heat flux can be calculated from [ab initio calculations](https://doi.org/10.1103/PhysRevLett.118.175901), [force fields](https://docs.lammps.org/compute_heat_flux.html), or [machine learning potentials](https://doi.org/10.1103/PhysRevB.108.L100302). 

## Example: LJ-Argon at 60 K

Here we use Lennard-Jones force fields as an example to calculate the thermal conductivity of Argon, for rapid prototyping the *ab initio* Green-Kubo workflow, including the noise reduction and extrapolation. Example files can be found in the [FHI-vibes Gitlab repository](https://gitlab.com/vibes-developers/vibes) at `example/green_kubo`.

The thermal conductivity are calculated from *NVE* trajectory. First we will pick a geometry from the *NVT* trajectory as the starting geometry for this *NVE* trajectory. 
Using command `vibes utils trajectory pick-sample -n 4500 trajectory.nc`, we pick the snapshot at step $4500$ from `trajectory.nc`. 
Together with `geometry.in.primitive` and `geometry.in.supercell`, we copy the `geometry.in.04500` to a new directory for preparing a *NVE* MD calculation. Here is an example of `md.in`,
```
[files]
geometry:                      geometry.in
primitive:                     geometry.in.primitive
supercell:                     geometry.in.supercell

[calculator]
name:                          lj

[calculator.parameters]
sigma:                         3.405
epsilon:                       0.010325
rc:                            13.0

[md]
driver:                        VelocityVerlet
timestep:                      4
maxsteps:                      20000
compute_stresses:              2
workdir:                       md

[md.kwargs]
logfile:                       md.log
```

To accelerate the calculation, we calculate the atomic stresses every $2$ time step.
After finishing the computation using command `vibes run md`, we use `vibes output md md/trajectory` to generate the trajectory file `trajectory.nc`.
Then we calculate the Green-Kubo thermal conductivity using command `vibes output gk --plot trajectory.nc`, a filter window size will be automatically chosen based on the lowest vibrational frequency from the MD simulation and the output shows like this:
```
Run aiGK output workflows for trajectory.nc
[GreenKubo]    Compute Prefactor:
[GreenKubo]    .. Volume:           4257.65  AA^3
[GreenKubo]    .. Temperature:        55.62  K
[GreenKubo]    -> Prefactor:      2.559e+10  W/mK / (eV/AA^2/fs)
[GreenKubo]    Estimate filter window size
[GreenKubo]    .. lowest vibrational frequency: 0.6500 THz
[GreenKubo]    .. corresponding window size:    1538.3590 fs
[GreenKubo]    .. window multiplicator used:    1.0000 fs
[filter]       Apply Savitzky-Golay filter with {'window_length': 385, 'polyorder': 1}
[GreenKubo]    Cutoff times (fs):
[GreenKubo]    [[1656.    0.    0.]
[GreenKubo]     [1328. 7592.    0.]
[GreenKubo]     [   0.    0. 5404.]]
[GreenKubo]    Kappa is:       0.253 +/- 0.066 W/mK
[GreenKubo]    Kappa^ab is: 
[GreenKubo]    [[ 1.192e-01 -2.103e-05 -1.030e-04]
[GreenKubo]     [ 2.268e-02  2.405e-01 -5.958e-05]
[GreenKubo]     [-1.109e-04 -7.210e-05  3.991e-01]]
.. write to greenkubo.nc
..    green kubo summary plotted to greenkubo_summary.pdf
```

??? info "greenkubo_summary.pdf"
    ![image](assets/greenkubo_summary.png)

The upper plot shows the heat flux auto-correlation function, with the thick line 
showing the curve after removing the non-contributing component and a filter. The 
lower plot shows the integration of the HFACF, thermal conductivity, of each Cartesian 
component and its diagonal mean and the standard error.

$$\kappa_{\rm mean} = 1/3 (\kappa_{xx} + \kappa_{yy} + \kappa_{zz})$$

Argon is cubic system showing space group $\rm F\bar{4}3m$, following Neumann’s principle, the 
symmetry operation reduces the thermal conductivity tensor to one independent coefficients 
with $\kappa_{xx} = \kappa_{yy} = \kappa_{zz}$. The error of the three components shown in 
the plot is arised from thermal fluctuation and statistic error.

### Finite-size correction

Green-Kubo theory calculate the thermal conductivity from a finite supercell, thus, the phonon wavelength greater than the cell size is missing in our calculation. Therefore, we need a finite-size correction for includes those long wavelength phonon modes, mostly acoustic modes.
The theory of the finite-size correction is introduced in [Carbogno PRL 2017](https://doi.org/10.1103/PhysRevLett.118.175901) and [Knoop PRB 2023](https://doi.org/10.1103/PhysRevB.107.224304), and the method is implemented in `FHI-vibes`. Here we will introduce the workflow for this feature.


#### Calculate the harmonic heat flux and finite-size extrapolation

Then, we calculate the harmonic heat flux using the force constant that we just prepared, and update the `trajectory.nc` with the following command:
```
vibes output md -fc output/FORCE_CONSTANTS --force md/trajectory.son
```

The finite-size extrapolation can be easily calculated with the following command:
```
vibes output gk --interpolate --plot trajectory.nc
```

You're expected to get something looks like this, which means you have done the right finite-size extrapolation.

```
[lattice_points] .. matched 108 positions in supercell and primitive cell in 0.470s
[symmetry]     reduce q-grid w/ 108 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  108/108
[symmetry]     .. q-points reduced from 108 to 20 points. in 0.079s
[force_constants] remap force constants
[force_constants] .. time elapsed: 0.122s
[force_constants] -> Symmetrize force constants.
[dynamical_matrix] Setup complete, eigensolution is unitary.
[gk.harmonic]  Get lifetimes by fitting to exponential
[gk.harmonic]  ** acf drops fast for s, q: 0, 0 set tau_sq = np.nan
[gk.harmonic]  ** acf drops fast for s, q: 1, 0 set tau_sq = np.nan
[gk.harmonic]  ** acf drops fast for s, q: 2, 0 set tau_sq = np.nan
[gk.harmonic]  .. time elapsed: 0.360s
[gk.interpolation] Interpolated l_sq at Gamma: [1.52 2.49 4.26]
[symmetry]     reduce q-grid w/ 64 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  64/64
[symmetry]     .. q-points reduced from 64 to 16 points. in 0.059s
[gk.interpolation]   4, Nq_eff =   0.59, kappa = 0.135 W/mK
[symmetry]     reduce q-grid w/ 216 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  216/216
[symmetry]     .. q-points reduced from 216 to 40 points. in 0.062s
[gk.interpolation]   6, Nq_eff =   2.00, kappa = 0.172 W/mK
[symmetry]     reduce q-grid w/ 512 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  512/512
[symmetry]     .. q-points reduced from 512 to 80 points. in 0.084s
[gk.interpolation]   8, Nq_eff =   4.74, kappa = 0.192 W/mK
[symmetry]     reduce q-grid w/ 1000 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  1000/1000
[symmetry]     .. q-points reduced from 1000 to 140 points. in 0.120s
[gk.interpolation]  10, Nq_eff =   9.26, kappa = 0.203 W/mK
[symmetry]     reduce q-grid w/ 1728 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  1728/1728
[symmetry]     .. q-points reduced from 1728 to 224 points. in 0.208s
[gk.interpolation]  12, Nq_eff =  16.00, kappa = 0.209 W/mK
[symmetry]     reduce q-grid w/ 2744 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  2744/2744
[symmetry]     .. q-points reduced from 2744 to 336 points. in 0.382s
[gk.interpolation]  14, Nq_eff =  25.41, kappa = 0.215 W/mK
[symmetry]     reduce q-grid w/ 4096 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  4096/4096
[symmetry]     .. q-points reduced from 4096 to 480 points. in 0.643s
[gk.interpolation]  16, Nq_eff =  37.93, kappa = 0.219 W/mK
[symmetry]     reduce q-grid w/ 5832 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  5832/5832
[symmetry]     .. q-points reduced from 5832 to 660 points. in 1.676s
[gk.interpolation]  18, Nq_eff =  54.00, kappa = 0.221 W/mK
[symmetry]     reduce q-grid w/ 8000 points
[symmetry]        |||||||||||||||||||||||||||||||||||||  8000/8000
[symmetry]     .. q-points reduced from 8000 to 880 points. in 1.858s
[gk.interpolation]  20, Nq_eff =  74.07, kappa = 0.225 W/mK
[gk.interpolation] Initial harmonic kappa value:       0.149 W/mK
[gk.interpolation] Fit intercept:                      0.247 W/mK
[gk.interpolation] Fit intercept - initial value:      0.098 +/- 0.001  W/mK
[gk.interpolation] Interpolated harm. kappa:           0.247 +/- 0.001 W/mK
[gk.interpolation] Correction^ab: 
[gk.interpolation] [[ 0.098 -0.003 -0.003]
[gk.interpolation]  [-0.003  0.098 -0.003]
[gk.interpolation]  [-0.003 -0.003  0.098]]
[gk.interpolation] Correction:                         0.098 +/- 0.001 W/mK
[gk.interpolation] Correction factor:                  1.655 +/- 0.008
[GreenKubo]    END RESULT: Finite-size corrected thermal conductivity
[GreenKubo]    Corrected kappa is:       0.350 +/- 0.066 W/mK
[GreenKubo]    Corrected kappa^ab (W/mK) is: 
[GreenKubo]    [[ 0.217 -0.003 -0.003]
[GreenKubo]     [ 0.019  0.338 -0.003]
[GreenKubo]     [-0.003 -0.003  0.497]]
.. write to greenkubo.nc
..    green kubo summary plotted to greenkubo_summary.png
.. interpolation summary plotted to greenkubo_interpolation.png
.. interpolation summary plotted to greenkubo_interpolation_fit.png
..      lifetime summary plotted to greenkubo_interpolation_lifetimes.png
```

and the following plots, helping to understand the extrapolation method.

??? info "greenkubo_interpolation_lifetimes.pdf"
    ![image](assets/greenkubo_interpolation_lifetimes.png)

This plot shows the phonon mode energy auto-correlation function. And the
analytical exponential decay correlation function.


??? info "greenkubo_interpolation.pdf"
    ![image](assets/greenkubo_interpolation.png)

The second plot show the thermal conductivity calculated from the Green-Kubo formula 
(as we showed before) and the BTE-like model at the commensurate q points, the phonon-phonon 
interaction is calculated from the real time auto-correlation function as shown in the 
left panel in the above plot.

??? info "greenkubo_interpolation_fit.pdf"
    ![image](assets/greenkubo_interpolation_fit.png)

The plots shows the interpolation and  extrapolation results. The thermal 
conductivity from harmonic model in recirpocal space is shown by the green square,
the x axis show $1/n_{q}$ with $n_{q}$ being the number of commensurate q points.
The blue dots shows the thermal conductivity results from the interpolated q grid.
As indicated by the orange line, the thermal conductivity increases with denser q grid.
By assuming the linear dependent of $\kappa_{n_q} \propto 1/n_q$, we can obtain the 
thermal conductivity at infinite q grid density $\kappa_{\rm bm-bulk}$ as the intercept 
of the linear fit. The finite-size correction $\Delta \kappa$ is calculated by,

$$\Delta\kappa = \kappa_{\rm hm-bulk} - \kappa_{\rm hm}$$




#### TODO
Adjusting the simulation time (for example double or half the simulation time) to test if the thermal conductivity increase.
Pick several snapshots from an NVT trajectory as the starting points as starting geometry, and run NVE trajectory from each starting points. Test the ensemble uncertainty and convergence.

