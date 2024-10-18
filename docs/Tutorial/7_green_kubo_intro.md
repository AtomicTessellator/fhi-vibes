# Green-Kubo theory

## *Ab initio* Green-Kubo theory

Green-Kubo theory is based on linear response theory, it calculates the thermal 
conductivity from equilibrium molecular dynamics. Green-Kubo theory relates the 
thermal conductivity from the integration of the real time autocorrelation function 
of heat flux $\boldsymbol{J}$,

$$\boldsymbol{\kappa}^{\alpha \beta}(T) = \frac{V}{k_B T^2} \int_{0}^{t_0} \langle \boldsymbol{J}^{\alpha}(t) \cdot \boldsymbol{J}^{\beta}(0) \rangle_T \ dt$$

where $\boldsymbol\kappa^{\alpha \beta}$ is the thermal conductivity tensor element 
in Cartesian coordinates of $\alpha, \beta$, with the volume of the supercell volume 
$V$, Boltzmann constant $k_B$ and temperature $T$.

In the density functional theory, the ab initio heat flux $\boldsymbol{J}_v$ 
(discarding the negligible convective term in solids) can be uniquely expressed from the 
virial theorm on the atomic contribution of stress tensor $\boldsymbol{\sigma}_{I}$
[[Carbogno2017](reference#Carbogno2017)]

$$\boldsymbol{J}_v(t) = \frac{1}{V}\sum_{I}\boldsymbol{\sigma}_{I}(t)\boldsymbol{\dot{R}}_I(t)$$

The atomic contribution to the stress tenor $\boldsymbol{\sigma}_I(t)$ and the virial heat flux 
$\boldsymbol{J}_v(t)$ will be calculated along the ab initio trajectory accounts for the full 
anharmonicity.


## Cutoff time and noise filtering

Limited by the supercell size and simulation time, the time correlation function suffer
heavy noise at its tail. Therefore, a cutoff time $t_0$ is chosen to avoid the 
statistical fluctuations. 

For a robust identification of the cutoff time, we first 
reduce noise by discarding the non-contributing term in atomic stress tensor
$\delta \boldsymbol{\sigma}_I(t) = \boldsymbol{\sigma}_I(t) - \langle \boldsymbol{\sigma}_I(t) \rangle$, 
and in heat flux $\delta \boldsymbol{J}(t) = \boldsymbol{J}(t) - \langle \boldsymbol{J}(t) \rangle$.

Then we apply a filter to remove the noise term that do not contribute to the 
integration function. The filtering window is automatically chosen based on the slowest
significant frequency $t_{\rm window} = 1/\omega_{\rm min}$, which is chosen to be the 
first peak in the vibrational density of states (VDOS).

Finally, the cutoff time $t_0$ is chosen from the first dip method, where the correlation function 
drops to zero.

## Finite-size correction

Due to the limited supercell size used in *ai*MD simulation, the vibration with longer 
wavelength than the supercell dimension is not included. This is particularly important 
for some harmonic system where the long-wavelength phonon mode contribute significantly 
to the heat transport.

### Harmonic mapping

To map the real space dynamics to the phonon picture which allows for interpolating in 
reciprocal space. We begin with the dynamical matrix in reciprocal space,

$$D_{IJ}^{\alpha\beta}(\boldsymbol{q}) = \sum_{L} \frac{1}{\sqrt{M_I M_J}} \Phi_{IJ}^{\alpha \beta} \exp(i\boldsymbol{q}(\boldsymbol{R}_I - \boldsymbol{R}_J - \boldsymbol{R}_L))$$

with $\Phi_{IJ}^{\alpha \beta}$ being the force constant of atom $I,J$ in Cartesian 
coordinates of $\alpha, \beta$, $M_I$ being the atomic mass and $\boldsymbol{R}_I$ 
being the reference position in the unit cell. $\boldsymbol{R}_L$ is a Bravais 
lattice vector. $\boldsymbol{q}$ is the commensurate wave vector. 

Note that in our *ab initio* molecular dynamics, the classical equation of motion is 
used for updating the position and momentum of nuclei, neglecting the quantum nuclei effect.
So we also will use classical equation of motion for the lattice dynamics, which is 
different from the quantum case. The solve of the dynamical matrix in classical case
can be written as an eigenvalue problem,

$$\sum_{J\beta}D^{\alpha\beta}_{IJ}(\boldsymbol{q}) \boldsymbol{e}_{s\boldsymbol{q}J} = \omega_{s\boldsymbol{q}}^2 \boldsymbol{e}_{s\boldsymbol{q}I}$$

which yields real eigenvalues $\omega_{s\boldsymbol{q}}^2$ and complex eigenvectors 
$\boldsymbol{e}_{s\boldsymbol{q}I}$. Using the eigenvector, we can map the position, 
$\boldsymbol{r}_I(t)$, and velocity, $\dot{\boldsymbol{r}}_I(t)$, in molecular dynamics 
to the normal coordinates $u_{s\boldsymbol{q}}(t)$ and momentum $p_{s\boldsymbol{q}}(t)$,

$$\begin{split} r_{sq}(t) = & \sum_{I}\sqrt{M_I} \cdot \boldsymbol{e}_{sqI} \cdot \boldsymbol{r}_I(t) \\ p_{sq}(t) = & \sum_{I} \frac{1}{\sqrt{M_I}} \cdot \boldsymbol{e}_{sqI} \cdot \dot{\boldsymbol{r}}_I(t) \end{split}$$

From here we can calculate time-dependent mode amplitude $a_{s\boldsymbol{q}}(t)$ and mode energy $E_{s\boldsymbol{q}}(t)$,

$$ a_{sq}(t) = \frac{1}{\sqrt{2}} \left( r_{sq}(t) + \frac{i}{\omega_{sq}} p_{sq}(t) \right)$$

$$ E_{sq}(t) = \omega^2_{sq}a^{\dagger}_{sq}(t)a_{sq}(t)$$

Using the mode-resolved energy, the harmonic energy flux $\boldsymbol{J}_{\rm{ha-q}}(t)$ can be defined as

$$\boldsymbol{J}_{\rm{ha-q}}(t) = \frac{1}{V} \sum_{sq} E_{sq}(t) \boldsymbol{v}_{sq}$$

For a classical harmonic system, we can use that $\langle E_{s\boldsymbol{q}}^2 \rangle = (k_B T)^2$ 
and neglecting the cross correlation between different modes. So the harmonic thermal 
conductivity is,

$$\boldsymbol{\kappa}^{\alpha \beta}_{\rm ha-\boldsymbol{q}}(T) = \frac{k_B}{V} \sum_{s\boldsymbol{q}} \boldsymbol{v}_{sq}^{\alpha} \boldsymbol{v}_{sq}^{\beta} \int_{0}^{t_0} \frac{\langle E_{s\boldsymbol{q}}(t) \cdot E_{s\boldsymbol{q}}(0) \rangle_T}{\langle E_{s\boldsymbol{q}}^2 \rangle} \ dt$$

By defining the phonon lifetime as the integral of the normalized time-dependent phonon energy,

$$\tau_{s\boldsymbol{q}} =  \int_{0}^{t_0} \frac{\langle E_{s\boldsymbol{q}}(t) \cdot E_{s\boldsymbol{q}}(0) \rangle_T}{\langle E_{s\boldsymbol{q}}^2 \rangle} \ dt$$

the above equation can be simplified to,

$$\boldsymbol{\kappa}^{\alpha \beta}_{\rm ha-\boldsymbol{q}}(T) = \frac{k_B}{V} \sum_{s\boldsymbol{q}} \boldsymbol{v}_{s\boldsymbol{q}}^{\alpha} \boldsymbol{v}_{s\boldsymbol{q}}^{\beta} \tau_{s\boldsymbol{q}}$$

### Lifetime interpolation

So far, we have calculated the harmonic thermal conductivity in commensurate 
$\boldsymbol{q}$ points in the supercell. This equation also allows us to interpolate 
the harmonic thermal conductivity to denser $\tilde{\boldsymbol{q}}$ grid. The phonon lifetime $\tau_{s}(\boldsymbol{q})$ scaling with 
$\omega_{s}^{-2}(\boldsymbol{q})$ which is rooted in basic phonon theory [[Pomeranchuk1941](references.md#Pomeranchuk1941),
[Herring1954](references.md#Herring1954)]. The frequency scaling can be written as,

$$\tau_s(\boldsymbol{q}) = \lambda_s(\boldsymbol{q}) \omega_{s}^{-2}(\boldsymbol{q})$$

The scaling factor is only weaky dependent on $\boldsymbol{q}$ vectors, so we linearly 
interpolate the $\lambda_s(\tilde{\boldsymbol{q}})$ as a function of $\tilde{\boldsymbol{q}}$ 
respectively for each mode $s$ in reciprocal space. 
The phonon frequency $\omega_{s}(\tilde{\boldsymbol{q}})$ and phonon group velocity
$\boldsymbol{v}_{s}(\tilde{\boldsymbol{q}})$ at dense $\boldsymbol{q}$ grid can be calculated
by Fourier interpolation in the dynamical matrix. 

For a given denser grid, the harmonic thermal conductivity can be obtained by

$$\boldsymbol{\kappa}^{\alpha \beta}_{\rm ha-int}(N_{\tilde{\boldsymbol{q}}}) = \frac{k_B}{V} \frac{N_{\boldsymbol{q}}}{N_{\tilde{\boldsymbol{q}}}} \sum_{s\tilde{\boldsymbol{q}}} \boldsymbol{v}_{s}^{\alpha}(\tilde{\boldsymbol{q}}) \boldsymbol{v}_{s}^{\beta}(\tilde{\boldsymbol{q}}) \tau_{s}(\tilde{\boldsymbol{q}})$$

### Thermal conductivity extrapolation

The bulk limit is reached at infinite grid density in reciprocal space. To access the 
value, the interpolated value is computed for an increasing density of $\boldsymbol{q}$ 
points. The interpolated value is approximately linear scaling with respective to 
the number of $\boldsymbol{q}$ points in each direction, $n_{\tilde{\boldsymbol{q}}} \equiv N_{\tilde{\boldsymbol{q}}}^{-1/3}$. So the bulk limit $\boldsymbol{\kappa}_{\rm ha-bulk}$ is 
obtained as the interception of the linear fitting.

Finally, the finite-size correction is calculated by the difference between the
harmonic thermal conductivity at the bulk limit and at the commensurate $\boldsymbol{q}$ points 
in the supercell.

$$\Delta \kappa = \kappa_{\rm ha-bulk} - \kappa_{ha-\boldsymbol{q}}$$



